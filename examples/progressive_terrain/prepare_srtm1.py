"""Download/extract/convert SRTM1 data from USGS:

https://www.usgs.gov/centers/eros/science/usgs-eros-archive-digital-elevation-shuttle-radar-topography-mission-srtm-1-arc?qt-science_center_objects=0#qt-science_center_objects
"""
import numpy as np
import glob
import cv2
import os
import tqdm
import wget
import zipfile
import argparse


def hgt_to_png(path):
    #n = 1201  # Change this to 3601 for SRTM1
    n = 3601
    with open(path, 'rb') as hgt_data:
        hgt = np.fromfile(hgt_data, np.dtype('>i2'), n * n).reshape((n, n)) / 10
        cv2.imwrite(path.replace('hgt', 'png'), hgt)


def download_srtm1_region01(path, min_lat=38, max_lat=49, min_lon=112, max_lon=124):
    os.makedirs(path, exist_ok=True)
    n_chunks = ((max_lon - min_lon) * (max_lat - min_lat))
    with tqdm.tqdm(total=n_chunks, desc='Downloading/extracting SRTM1...') as pbar:
        for i in range(min_lat, max_lat + 1):
            for j in range(min_lon, max_lon + 1):
                name = f'N{i}W{j}.hgt.zip'
                url = f'https://dds.cr.usgs.gov/srtm/version2_1/SRTM1/Region_01/{name}'
                expected = os.path.join(path, name)
                if not os.path.exists(expected):
                    try:
                        filename = wget.download(url, out=path)
                        zipped = zipfile.ZipFile(filename)
                        zipped.extractall(path)
                        zipped.close()
                        hgt = os.path.join(path, name.replace('.zip', ''))
                        hgt_to_png(hgt)
                    except:
                        print(f'Failed to prepare: {name}')
                pbar.update(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Make PNGs from SRTM1 data')
    parser.add_argument('--output', type=str, default='./srtm1/',
                        help='Where to put all data')
    parser.add_argument('--min_lat', type=int, default=38,
                        help='Minimum latitude to prepare')
    parser.add_argument('--max_lat', type=int, default=38,
                        help='Maximum latitude to prepare')
    parser.add_argument('--min_lon', type=int, default=112,
                        help='Minimum longitude to prepare')
    parser.add_argument('--max_lon', type=int, default=112,
                        help='Maximum longitude to prepare')
    args = parser.parse_args()

    download_srtm1_region01(args.output,
                            args.min_lat, args.max_lat,
                            args.min_lon, args.max_lon)
