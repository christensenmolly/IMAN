# DESCRIPTION:
# Скрипт для оценки стандартных отклонений параметров galfit-модели с помощью Монте-Карло моделирования фона неба
# (фон неба считается константным, его уровеньизвлекается из нормального распределения с нулевым средним и
# заданным стандартным отклонением, в качестве стандартного отелонения предполагается брать rms фона неба).
# USAGE: python estimate_galfit_residuals.py [galfit_file_absolute_path] --sky_rms [sky_rms] --n [n]
# --dist_modulus [dist_modulus] --sec_in_pixel [sec_in_pixel] --kpc_in_second [kpc_in_second]

import argparse
import numpy as np
import os
import pandas as pd
import shutil
import subprocess

from os.path import basename, join, exists, isdir
from statistics import stdev

SKY_CONFIG_TEMPLATE = ("# Component: sky \n", "0) sky #  Component type\n",
                       "1) {}    0 #  Sky background at center of fitting region [ADUs]\n",
                       "2) 0.000e+00  0 #  dsky/dx (sky gradient in x)     [ADUs/pix]\n",
                       "3) 0.000e+00  0 #  dsky/dy (sky gradient in y)     [ADUs/pix]\n",
                       "Z) 0            #  Skip this model in output image?  (yes=1, no=0)\n")
STOP_LINE = "# INITIAL FITTING PARAMETERS"
REQUIRED_PARAMETERS = ["A)", "D)", "F)"]
KEY_STRING = "0) sky"


def copy_galfit_required(galfit_inp: str, work_directory: str):
    """
    Считывает из galfit.inp необходимые файлы и копирует их в work_directory.
    Необходимыми являются: A) Input data image (FITS file),
                           D) Input PSF image and (optional) diffusion kernel,
                           F) Bad pixel mask (FITS image or ASCII coord list).
    :param galfit_inp: абсолютный путь до входного galfit-файла,
           !!! предполагается, что вышеперечисленные необходимые файлы
           содержаться в той же директории, что и galfit.inp
    :param work_directory: абсолютный путь до рабочей директории,
    """
    with open(galfit_inp, "r", encoding="utf8") as handler:
        galfit_lines = [s.strip() for s in handler.readlines()]
    for line in galfit_lines:
        if line == STOP_LINE:
            break
        if line and line[0] != "#" and line.split()[0] in REQUIRED_PARAMETERS:
            required_name = line.split()[1]
            required_path = join(os.path.dirname(galfit_inp), required_name)
            required_dst_path = join(work_directory, required_name)
            shutil.copy(required_path, required_dst_path)


def modify_galfit_inp(galfit_inp: str, sky_rms: float):
    """
    Модифицирует копию входного galfit-файла в рабочей директории,
    заменяя константный уровень неба на уровень из
    центрированного нормального распределения с заданным стандартным отклонением
    :param galfit_inp: абсолютный путь до копии входного galfit-файла
           в рабочей директории, который предполагается изменить
    :param sky_rms: rms фона неба, который берётся в качестве стандартного отклонения
           нормального распределения
    """
    with open(galfit_inp, "r", encoding="utf8") as reader:
        galfit_lines = reader.readlines()

    i = 0
    while (i + 1) < len(galfit_lines) and KEY_STRING not in galfit_lines[i]:
        i += 1
    if (i + 1) < len(galfit_lines):
        words = galfit_lines[i + 1].split()
        words[1] = str(np.random.normal(scale=sky_rms))
        galfit_lines[i + 1] = " ".join(words) + "\n"
    else:
        sky_config = list(SKY_CONFIG_TEMPLATE)
        sky_config[1] = sky_config[1].format(str(np.random.normal(scale=sky_rms)))
        galfit_lines.extend(sky_config)

    with open(galfit_inp, "w", encoding="utf8") as writer:
        writer.writelines(galfit_lines)


def run_galfit(galfit_inp: str, sky_rms: float, n=1, show_galfit_output=False):
    """
    Запускает galfit-декомпозицию n раз с разными уровнями фона неба
    :param galfit_inp: абсолютный путь до копии входного galfit-файла,
    :param sky_rms: rms фона неба, который берётся в качестве стандартного отклонения
           нормального распределения
    :param n: количество запусков galfit-декомпозиции
    :param show_galfit_output: булевый флаг, если True, то выводится выходной поток
           galfit-декомпозиции
    """
    galfit_inp_name = basename(galfit_inp)

    np.random.seed(42)
    for i in range(n):
        print("Decomposition number {} starts.".format(i))
        modify_galfit_inp(galfit_inp, sky_rms)
        process = subprocess.run(["galfit", galfit_inp_name], check=True,
                                 stdout=subprocess.PIPE, universal_newlines=True)
        if show_galfit_output:
            print(process.stdout)


def get_galfit_number(n: int) -> str:
    """
    Создаёт имя выходного galfit-файла по его номеру (n),
    предполагается, что 1 <= n <= 99.
    :param n: номер выходного galfit-файла
    :return: имя выходного galfit-файла
    """
    if n < 10:
        return "galfit.0{}".format(n)
    else:
        return "galfit." + str(n)


def get_data(n: int) -> dict:
    """
    Из полученных в ходе n декомпозиций выходных galfit-файлов извлекает
    значения параметров и помещает их в словарь, в котором ключи -- номера строк в
    galfit-файлах, значения -- списки значений параметров;
    предполагается, что 1 <= n <= 99.
    :param n: количество проведённых galfit-декомпозиций
    :return: словарь c значениями параметров
    """
    # Создаём словарь значений, где ключ -- номер строки,
    # значение -- список пустых списков
    # (длина списка соответствует количеству параметров)
    galfit_dict = dict()
    i = 1
    while i <= n and not exists(get_galfit_number(i)):
        i += 1
    if i <= n:
        galfit_oup = get_galfit_number(i)
    else:
        raise RuntimeError("No galfit output files have been created.")
    with open(galfit_oup, "r", encoding="utf8") as reader:
        galfit_lines = reader.readlines()
    i = 0
    while i < len(galfit_lines):
        if "# Component number:" in galfit_lines[i]:
            i += 1
            while galfit_lines[i].strip():
                elements = galfit_lines[i].partition("#")[0].split()
                if len(elements) > 2:
                    n_values = len(elements) // 2
                    galfit_dict[i] = [[] for i in range(n_values)]
                i += 1
        else:
            i += 1

    # Заполняем словарь соответствующими значениями
    for i in range(1, n + 1):
        galfit_oup = get_galfit_number(i)
        if exists(galfit_oup):
            with open(galfit_oup, "r", encoding="utf8") as reader:
                galfit_lines = reader.readlines()
            for key, value in galfit_dict.items():
                n_param = len(value)
                elements = galfit_lines[key].partition("#")[0].split()
                for j in range(1, n_param + 1):
                    value[j - 1].append(float(elements[j]))
        else:
            print("!!!Attention!!!\n{} has not been created.".format(galfit_oup))
    return galfit_dict


def process_line(galfit_line: str, max_len: int, value: list) -> str:
    """
    Формирует строку с параметрами в galfit.stat
    (итоговый файл, генерируемый скриптом) по самой исходной строке и
    значениям рассчитанных стандартных отклонений параметров,
    описываемых в этой строке, например:
    1) 414.7755 433.7625 0 0  #  Position x, y
    преобразуется в
    1) 414.7755 433.7625 0 0  #  Position x, y  0.000e+00  0.000e+00
    :param galfit_line: исходная galfit-строка
    :param max_len: максимальная длина строки в исходном galfit-файле
    :param value: список списков значений параметров
    :return: строка для galfit.stat
    """
    galfit_line = galfit_line.strip()
    length = len(galfit_line)
    result_line = galfit_line + " " * (max_len - length + 2)
    for elem in value:
        result_line += " {0:8.3e} ".format(stdev(elem))
    return result_line + "\n"


def write_statistics(galfit_inp: str, galfit_dict: dict):
    """
    По словарю значений параметров galfit-моделей создаёт файл
    galfit.stat с оценкой стандартных отклонений параметров модели.
    :param galfit_inp: абсолютный путь до копии входного galfit-файла
           в рабочей директории
    :param galfit_dict: словарь c значениями параметров
    """
    with open(galfit_inp, "r", encoding="utf8") as reader:
        galfit_lines = reader.readlines()

    max_len = max([len(line) for line in galfit_lines[20:]]) - 1
    for key, value in galfit_dict.items():
        galfit_lines[key] = process_line(galfit_lines[key], max_len, value)

    with open("galfit.stat", "w", encoding="utf8") as writer:
        writer.writelines(galfit_lines)


def write_tsv_statistics(galfit_inp: str, galfit_dict: dict, distance_modulus: float,
                         seconds_in_pixel: float, kpc_in_second: float):
    """
    По словарю значений параметров galfit-моделей создаёт tsv-таблицу:
    m	            M	    r_e[pix]	r_e[kpc]	  n
    $15 \pm 0.3$	$-22$	$6.6 \pm 5$	$6.7 \pm 5.5$	$6.5 \pm 4.2$
    ...
    galfit.stat с оценкой стандартных отклонений параметров модели.
    :param galfit_inp: абсолютный путь до копии входного galfit-файла
           в рабочей директории
    :param galfit_dict: словарь c значениями параметров
    :param distance_modulus: разность между наблюдаемой и абсолютной
           звёздными величинами
    :param seconds_in_pixel: количество угловых секунд в пикселе изображения
    :param kpc_in_second: количество kpc в угловой секунде
    """
    app_m = []
    abs_M = []
    r_e_pix = []
    r_e_kpc = []
    sersic_n = []
    with open(galfit_inp, "r", encoding="utf8") as reader:
        galfit_lines = [line.strip() for line in reader.readlines()]
        for key, value in galfit_dict.items():
            line = galfit_lines[key]
            if line.startswith("3)") and "Integrated magnitude" in line:
                app_m.append("${} \\pm {:.4f}$".format(line.split()[1],
                                                       round(stdev(value[0]), 4)))
                abs_M.append("${:.4f}$".format(round(float(line.split()[1]) -
                                                     distance_modulus), 4))
            if line.startswith("4)"):
                if "effective radius" not in line:
                    r_e_pix.append(" ")
                    r_e_kpc.append(" ")
                else:
                    r_e_pix.append("${} \\pm {:.4f}$".format(line.split()[1],
                                                             round(stdev(value[0]), 4)))
                    dev = stdev(value[0]) * seconds_in_pixel * kpc_in_second
                    val = float(line.split()[1]) * seconds_in_pixel * kpc_in_second
                    r_e_kpc.append("${:.4f} \\pm {:.4f}$".format(round(val, 4),
                                                                 round(dev, 4)))
            if line.startswith("5)"):
                if "Sersic index" not in line:
                    sersic_n.append(" ")
                else:
                    sersic_n.append("${} \\pm {:.4f}$".format(line.split()[1],
                                                              round(stdev(value[0]), 4)))
    stat_df = pd.DataFrame(list(zip(app_m, abs_M, r_e_pix, r_e_kpc, sersic_n)),
                           columns=["m", "M", "r_e[pix]", "r_e[kpc]", "n"])
    stat_df.to_csv("galfit_stat.tsv", sep="\t", encoding="utf8", index=False)


def main(galfit_inp: str, sky_rms: float, n: int, distance_modulus: float,
         seconds_in_pixel: float, kpc_in_second: float):
    """
    Главная управляющая функция модуля.
    По galfit-файлу вычисляет оценки стандартных отклонений параметров декомпозиции,
    используя Монте-Карло моделирование. Полученные оценки записываются в файл
    galfit.stat и таблицу galfit_stat.tsv
    :param galfit_inp: путь до исходного galfit-файла
    :param sky_rms: rms фона неба
    :param n: количсетво galfit-декомпозиций
    :param distance_modulus: разность между наблюдаемой и абсолютной
           звёздными величинами
    :param seconds_in_pixel: количество угловых секунд в пикселе изображения
    :param kpc_in_second: количество kpc в угловой секунде
    """
    if n < 2:
        raise ValueError("To estimate stddev we must have 2 or more samples (n > 1)")

    work_directory = join(os.path.dirname(galfit_inp), "residuals")
    if exists(work_directory):
        if isdir(work_directory):
            shutil.rmtree(work_directory)
        else:
            os.remove(work_directory)
    os.mkdir(work_directory)

    dist_galfit_path = join(work_directory, "galfit.inp")
    shutil.copy(galfit_inp, dist_galfit_path)
    copy_galfit_required(galfit_inp, work_directory)

    os.chdir(work_directory)
    run_galfit(dist_galfit_path, sky_rms, n)
    galfit_dict = get_data(n)
    write_statistics(dist_galfit_path, galfit_dict)
    write_tsv_statistics(dist_galfit_path, galfit_dict, distance_modulus,
                         seconds_in_pixel, kpc_in_second)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Monte Carlo modeling")
    parser.add_argument("galfit_input_file", help="Galfit input file absolute path")
    parser.add_argument("--sky_rms", nargs='?', help="Sky rms", type=float)
    parser.add_argument("--n", nargs='?', help="Number of galfit decompositon", type=int, default=2)
    parser.add_argument("--dist_modulus", nargs='?',
                        help="Difference between the apparent magnitude and the absolute magnitude", type=float)
    parser.add_argument("--sec_in_pixel", nargs='?', help="Second of arc in one pixel", type=float)
    parser.add_argument("--kpc_in_second", nargs='?', help="Kpc in second of arc", type=float)

    args = parser.parse_args()

    main(args.galfit_input_file, args.sky_rms, args.n, args.dist_modulus, args.sec_in_pixel, args.kpc_in_second)
