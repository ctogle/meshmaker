from noise import pnoise2, snoise2
import numpy as np
import cv2
from skimage import filters, segmentation
from skimage.measure import label, regionprops
from meshmaker.vec3 import vec3
from meshmaker.tform import tform


def steepness(img):
    xgrad, ygrad = np.gradient(img)
    grad_mag = np.sqrt(xgrad ** 2 + ygrad ** 2)
    #flat_mask = grad_mag <= filters.threshold_otsu(grad_mag)
    return grad_mag


def normalize(img):
    if img.max() > 0:
        return (img - img.min()) / img.max()
    else:
        return img


def dilate(img, kernal=(8, 8), r=10):
    return cv2.dilate(img, kernal, r)


def segment_regions(img):
    clean_border = segmentation.clear_border(img)
    labeled = label(clean_border)
    regions = regionprops(labeled)
    regions = sorted(regions, key=(lambda r: r.area), reverse=True)
    masks = [(labeled == r.label) for r in regions]
    fps = []
    for m, r in zip(masks, regions):
        s = vec3(1.0 / m.shape[0], 1.0 / m.shape[1], 1)
        fp = [(s * p) for p in pixel_polygon(m)]
        fps.append(fp)
    return fps, masks, regions


def perlin(resolution, frequency, octaves, alpha):
    img = np.ndarray((resolution, resolution))
    for y in range(resolution):
        for x in range(resolution):
            z = snoise2(x / frequency, y / frequency, octaves)
            img[y, x] = ((1 + z) / 2) ** alpha
    for i in range(3):
        img = cv2.blur(img, (8, 8))
    return img


def proximal(resolution, loop, f):
    contour = np.array([[[int(p.y), int(p.x)]] for p in loop])
    img = np.zeros((resolution, resolution))
    for i in range(resolution):
        for j in range(resolution):
            d = cv2.pointPolygonTest(contour, (j, i), True)
            img[j, i] = f(d)
    return img


def pixel_polygon(mask):
    pixels = mask.astype(np.uint8) * 255
    pixels = dilate(pixels, (8, 8), 10)
    ret, thresh = cv2.threshold(pixels, 127, 255, 0)
    contours, hierarchy = cv2.findContours(thresh, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
    if contours:
        vs = cv2.approxPolyDP(contours[0], 0.1, True)
        ps = [vec3(x, y, 0) for x, y in vs[:, 0]]
        return ps
