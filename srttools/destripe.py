import numpy as np
try:
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

from .fit import mad


def mask_zeros(image, expo=None, npix_tol=None):
    """Mask the lines containing zeros in the image.

    Parameters
    ----------
    image : 2d array
        Input image
    npix_tol : int
        Number of tolerated pixels with value 0

    Returns
    -------
    masked_image : 2d array
        The masked image
    mask : 2d array
        The boolean mask to obtain masked_image from mask

    Examples
    --------
    >>> import numpy as np
    >>> img = [[0, 1, 1], [0, 1, 1], [1, 1, 1]]
    >>> masked_image, mask = mask_zeros(img, expo=img, npix_tol=1)
    >>> np.all(masked_image == [[1, 1], [1, 1], [1, 1]])
    True
    >>> np.all(mask == [[False, True, True], [False, True, True],
    ...                 [False, True, True]])
    True
    >>> masked_image, mask = mask_zeros(img, npix_tol=2)
    >>> np.all(masked_image == img)
    True
    >>> img = [[0, 0, 0], [1, 1, 1], [1, 1, 1]]
    >>> masked_image, mask = mask_zeros(img, npix_tol=1)
    >>> np.all(masked_image == [[1, 1, 1], [1, 1, 1]])
    True
    """
    image = np.asarray(image)
    mask = np.ones(image.shape, dtype=bool)
    if npix_tol is None:
        return image, mask

    if expo is None:
        expo = image
    expo = np.asarray(expo)
    good_hor = 0
    for i in range(expo.shape[0]):
        line = expo[i, :]
        if len(line[line == 0]) > npix_tol:
            mask[i, :] = False
        else:
            good_hor += 1

    good_ver = 0
    for i in range(expo.shape[1]):
        line = expo[:, i]
        if len(line[line == 0]) > npix_tol:
            mask[:, i] = False
        else:
            good_ver += 1

    masked_image = image[mask].reshape((good_hor, good_ver))
    return masked_image, mask


def clip_and_smooth(img, clip_sigma=3, smooth_window=10, direction=0):
    """
    Examples
    --------
    >>> img = np.zeros((2,2))
    >>> np.all(clip_and_smooth(img, smooth_window=(5, 5)) == img)
    True
    >>> img = np.array([[0, 0], [1, 1]])
    >>> np.all(clip_and_smooth(img, direction=0) == img)
    True
    >>> img = np.array([[0, 1], [0, 1]])
    >>> np.all(clip_and_smooth(img, direction=1) == img)
    True
    >>> img = np.array([[1, 1.], [8., 1]])
    >>> np.allclose(clip_and_smooth(img, clip_sigma=1, smooth_window=0),
    ...             [[1, 1], [3.0310889132455352, 1]])
    True
    """
    from scipy.ndimage import gaussian_filter, gaussian_filter1d
    import collections
    if img.shape[0] * img.shape[0] > 100:
        rms = mad(img.flatten())
    else:
        rms = np.std(img.flatten())

    median = np.median(img)
    bad = img - median > clip_sigma * rms
    img[bad] = clip_sigma * rms
    bad = median - img > clip_sigma * rms
    img[bad] = - clip_sigma * rms

    if smooth_window == 0:
        pass

    elif isinstance(smooth_window, collections.Iterable):
        img = gaussian_filter(img, np.array(smooth_window) / 5)
    else:
        img = gaussian_filter1d(img, smooth_window / 5,
                                axis=np.logical_not(direction))
    return img


def basket_weaving(img_hor, img_ver, clip_sigma=3, niter_max=10,
                   expo_hor=None, expo_ver=None, window_shape='hanning'):
    """Basket-Weaving algorithm from Mueller et al. 1707.05573v6."""
    it = 1
    if expo_hor is None:
        expo_hor = np.ones_like(img_hor)
    if expo_ver is None:
        expo_ver = np.ones_like(img_ver)
    img_hor = np.copy(img_hor)
    img_ver = np.copy(img_ver)
    width = np.max(img_hor.shape)

    while it <= niter_max:
        window = width // 2 - 4 * it
        if window < 4:
            break

        diff = img_hor - img_ver
        diff = clip_and_smooth(diff, clip_sigma=clip_sigma,
                               smooth_window=(0., window))

        img_hor = img_hor - diff

        diff = img_ver - img_hor
        diff = clip_and_smooth(diff, clip_sigma=clip_sigma,
                               smooth_window=(window, 0.))

        img_ver = img_ver - diff
        it += 1

    img_final = img_ver * expo_ver + img_hor * expo_hor
    expo = expo_hor + expo_ver

    good = expo > 0
    img_final[good] = img_final[good] / expo[good]
    return img_final


def destripe_wrapper(image_hor, image_ver, alg='basket-weaving',
                     niter=10, expo_hor=None, expo_ver=None,
                     npix_tol=None, clip_sigma=3, label="img"):
    if expo_hor is None or expo_ver is None:
        image_mean = (image_hor + image_ver) / 2
        expo_hor = expo_ver = np.ones_like(image_mean)
        masked_image, mask = mask_zeros(image_mean, npix_tol=npix_tol)
    else:
        image_mean = \
            (image_hor*expo_hor + image_ver*expo_ver) / (expo_hor + expo_ver)
        masked_image, mask = mask_zeros(image_mean, expo_hor + expo_ver,
                                        npix_tol=npix_tol)

    if HAS_MPL:
        fig = plt.figure()
        plt.imshow(image_hor[mask].reshape(masked_image.shape))
        plt.savefig(label + '_hor.png')
        plt.imshow(image_ver[mask].reshape(masked_image.shape))
        plt.savefig(label + '_ver.png')
        diff_img = image_ver[mask] - image_hor[mask]
        plt.imshow(diff_img.reshape(masked_image.shape))
        plt.savefig(label + '_diff.png')
        plt.close(fig)

        fig = plt.figure()
        plt.imshow(expo_hor[mask].reshape(masked_image.shape))
        plt.savefig(label + '_expoh.png')
        plt.imshow(expo_ver[mask].reshape(masked_image.shape))
        plt.savefig(label + '_expov.png')
        plt.imshow(image_mean[mask].reshape(masked_image.shape))
        plt.savefig(label + '_initial.png')
        plt.close(fig)

    image_mean[mask] = \
        basket_weaving(image_hor[mask].reshape(masked_image.shape),
                       image_ver[mask].reshape(masked_image.shape),
                       niter_max=niter,
                       expo_hor=expo_hor[mask].reshape(masked_image.shape),
                       expo_ver=expo_ver[mask].reshape(masked_image.shape),
                       clip_sigma=clip_sigma
                       ).flatten()

    if HAS_MPL:
        plt.imshow(image_mean[mask].reshape(masked_image.shape))
        plt.savefig(label + '_destr.png')

    if alg == 'basket-weaving':
        return image_mean
