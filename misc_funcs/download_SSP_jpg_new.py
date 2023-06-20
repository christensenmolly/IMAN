#!/bin/env python3

import argparse
import tarfile
import subprocess
import tempfile
import getpass
import os, os.path
import contextlib
import logging ; logging.basicConfig(level=logging.INFO)
try:
    import pyfits
except:
    import astropy.io.fits as pyfits
import numpy
import PIL.Image

TOP_PAGE = 'https://hsc-release.mtk.nao.ac.jp/das_cutout/pdr3/'
API = 'https://hsc-release.mtk.nao.ac.jp/das_cutout/pdr3/cgi-bin/cutout'

def main(ra, dec, radius, name, outDir, user, password, filters=['HSC-I', 'HSC-R', 'HSC-G'], rerun='any', color='hsc'):
    # radius is in arcmin
    
    checkPassword(user, password)

    #coords, outs = loadCoords(args.input)
    coords = []
    outs = []
    for k in range(len(ra)):
        coords.append([ra[k], dec[k]])
        outs.append(name[k])    
    
    mkdir_p(outDir)
    

    
    batchSize = 600
    for batchI, batchCoords in enumerate(batch(coords, batchSize)):
        with requestFileFor(batchCoords, filters, radius, rerun) as requestFile:
            tarMembers = queryTar(user, password, requestFile)
            for i, rgb in rgbBundle(tarMembers):
                j = batchI * batchSize + i
                outFile = os.path.join(outDir, outs[j])
                logging.info('-> {}'.format(outFile))
                makeColorPng(rgb, outFile, color)





def loadCoords(input):
    import re
    comment = re.compile('\s*(?:$|#)')
    num = 1
    coords = []
    outs = []
    for line in input:
        if comment.match(line):
            continue
        cols = line.split()
        if len(cols) == 2:
            ra, dec = cols
            out = '{}.png'.format(num)
        else:
            ra, dec, out = cols
        ra = float(ra)
        dec = float(dec)
        num += 1
        coords.append([ra, dec])
        outs.append(out)
    return coords, outs


@contextlib.contextmanager
def requestFileFor(coords, filters, radius, rerun):
    with tempfile.NamedTemporaryFile() as tmp:
        tmp.write('#? filter ra dec sw sh rerun\n'.encode('utf-8'))
        for k in range(len(coords)):
            coord = coords[k]
            fov = '%iasec' % (int(radius[k]*60.*2.)) # '30asec'
            for filterName in filters:
                tmp.write('{} {} {} {} {} {}\n'.format(filterName, coord[0], coord[1], fov, fov, rerun).encode())
        tmp.flush()
        yield tmp.name


def batch(arr, n):
    i = 0
    while i < len(arr):
        yield arr[i : i + n]
        i += n


def rgbBundle(files):
    mktmp = tempfile.NamedTemporaryFile
    with mktmp() as r, mktmp() as g, mktmp() as b:
        lastObjNum = 0
        rgb = {}
        for info, fileObj in files:
            resNum = int(os.path.basename(info.name).split('-')[0])
            objNum = (resNum - 2) // 3
            if lastObjNum != objNum:
                yield lastObjNum, rgb
                rgb.clear()
                lastObjNum = objNum
            ch = 'gbr'[resNum % 3]
            dst = locals()[ch]
            copyFileObj(fileObj, dst)
            rgb[ch] = dst.name
        yield objNum, rgb



def copyFileObj(src, dst):
    dst.seek(0)
    dst.write(src.read())
    dst.truncate()


def checkPassword(user, password):
    with tempfile.NamedTemporaryFile() as netrc:
        netrc.write('machine hsc-release.mtk.nao.ac.jp login {} password {}'.format(user, password).encode('ascii'))
        netrc.flush()
        httpCode = subprocess.check_output(['curl', '--netrc-file', netrc.name, '-o', os.devnull, '-w', '%{http_code}', '-s', TOP_PAGE]).strip()
        if httpCode == b'401':
            raise RuntimeError('Account or Password is not correct')


def queryTar(user, password, requestFile):
    with tempfile.NamedTemporaryFile() as netrc:
        netrc.write('machine hsc-release.mtk.nao.ac.jp login {} password {}'.format(user, password).encode('ascii'))
        netrc.flush()

        pipe = subprocess.Popen([
            'curl', '--netrc-file', netrc.name,
            '--form', 'list=@{}'.format(requestFile),
            '--silent',
            API,
        ], stdout=subprocess.PIPE)

        with tarfile.open(fileobj=pipe.stdout, mode='r|*') as tar:
            while True:
                info = tar.next()
                if info is None: break
                logging.info('extracting {}...'.format(info.name))
                f = tar.extractfile(info)
                yield info, f
                f.close()

        pipe.wait()


def makeColorPng(rgb, out, color):
    if len(rgb) == 0:
        return

    with pyfits.open(list(rgb.values())[0]) as hdul:
        template = hdul[1].data

    layers = [numpy.zeros_like(template) for i in range(3)]
    for i, ch in enumerate('rgb'):
        if ch in rgb:
            with pyfits.open(rgb[ch]) as hdul:
                x = scale(hdul[1].data, hdul[0].header['FLUXMAG0'])
                layers[i] = x

    if color == 'hsc':
        layers = hscColor(layers)
    elif color == 'sdss':
        layers = sdssColor(layers)
    else:
        assert False

    layers = numpy.array(layers)
    layers[layers < 0] = 0
    layers[layers > 1] = 1
    layers = layers.transpose((1, 2, 0))[::-1, :, :]

    layers = numpy.array(255 * layers, dtype=numpy.uint8)
    img = PIL.Image.fromarray(layers)
    img.save(out)


def scale(x, fluxMag0):
    mag0 = 19
    scale = 10 ** (0.4 * mag0) / fluxMag0
    x *= scale
    return x


def hscColor(rgb):
    u_min = -0.05
    u_max = 2. / 3.
    u_a = numpy.exp(10.)
    for i, x in enumerate(rgb):
        x = numpy.arcsinh(u_a*x) / numpy.arcsinh(u_a)
        x = (x - u_min) / (u_max - u_min)
        rgb[i] = x
    return rgb


def sdssColor(rgb):
    u_a = numpy.exp(10.)
    u_b = 0.05
    r, g, b = rgb
    I = (r + g + b) / 3.
    for i, x in enumerate(rgb):
        x = numpy.arcsinh(u_a * I) / numpy.arcsinh(u_a) / I * x
        x += u_b
        rgb[i] = x
    return rgb


def mkdir_p(d):
  if not os.path.exists(d):
    try:
        os.makedirs(d)
    except:
        pass


if __name__ == '__main__':
    main()
