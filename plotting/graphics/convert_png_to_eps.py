#!/usr/bin/python
import argparse
import subprocess
import os

def rewrite_to_eps(input_image, output_image=None):
    '''
    arr_img = plt.imread(input_image, format='png')

    plt.axis("off")
    plt.imshow(arr_img)
    plt.savefig(output_image, transparent = False, dpi=300, bbox_inches='tight', pad_inches=0.0)    
    '''
    if output_image is None:
        output_image = input_image.split('.')[0] + '.eps'
    
    subprocess.call('/usr/bin/convert %s eps2:%s' % (input_image, output_image), shell=True)
    subprocess.call('epstopdf %s' % (output_image), shell=True)
    subprocess.call('pdftops -eps %s' % (output_image.split('.eps')[0] + '.pdf'), shell=True)
    os.remove(output_image.split('.eps')[0] + '.pdf')

    #subprocess.call('/usr/bin/convert %s eps3:%s' % (output_image, output_image), shell=True)
    
    #subprocess.call('/usr/bin/convert %s %s' % (input_image, output_image.split('.eps')[0] + '.pdf'), shell=True)
    #subprocess.call('pdftops -eps %s' % (output_image.split('.eps')[0] + '.pdf'), shell=True)
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Masking")
    parser.add_argument("input_image", help="Input image")
    parser.add_argument("output_image", nargs='?', const=1, help="Optional: Ouput image name",type=str, default=None) 
    args = parser.parse_args()

    input_image = args.input_image
    output_image = args.output_image   

    rewrite_to_eps(input_image, output_image)