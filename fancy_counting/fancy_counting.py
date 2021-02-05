import numpy as np
from scipy import special


def cdf_gauss(x):
    return (1 + special.erf(x)) / 2


# probabilities is 1D array (or list) with the props for a single his or miss (per bin)
# set value 0 for start
# returns <x>, <x**2>
# O(2**n) worst, if no 1s or 0s
def calc_mean_var(probabilities, value=0, k=0, prob=None):
    # if prob is None initilize
    if prob is None:
        left = calc_mean_var(probabilities, value=1, k=1, prob=probabilities[k])
        right = calc_mean_var(probabilities, value=0, k=1, prob=1-probabilities[k])
        return left[0] + right[0], left[1] + right[1]
    if prob == 0:
        return 0, 0
    if k == len(probabilities):
        return value*prob, value**2*prob
    else:
        left = calc_mean_var(probabilities, value=value + 1, k=k + 1, prob=prob * probabilities[k])
        right = calc_mean_var(probabilities, value=value, k=k + 1, prob=prob * (1-probabilities[k]))
        return left[0] + right[0], left[1] + right[1]


# O(n**2)
def calc_mean_var_fast(probabilities):
    cum_probs = [1]
    sure_thing = 0
    for p in probabilities:
        if p == 0:
            pass
        elif p == 1:
            # insert is O(n) !
            # cum_probs.insert(0, 0)
            # all things that will definitely happen, are summed up in one int, that will be added to x later
            sure_thing += 1
        else:
            # see block comment below
            last = cum_probs[-1]
            cum_probs = [cum_probs[0] * (1 - p) if k == 0 else cum_probs[k] * (1 - p) + cum_probs[k - 1] * p for k in range(len(cum_probs))]
            cum_probs.append(last * p)
    # mean = sum([i * cum_probs[i] for i in range(len(cum_probs))])
    # e_x_sqared = sum([i**2 * cum_probs[i] for i in range(len(cum_probs))])
    # returns <x> = sum(x_i * p(x_i)), <x**2> = sum(x_i**2 * p(x_i))
    return sum([(i + sure_thing) * cum_probs[i] for i in range(len(cum_probs))]), sum([(i + sure_thing)**2 * cum_probs[i] for i in range(len(cum_probs))])


# n-bins = len(bins)-1, because upper and lowest boundary are given
def build_prob_matrix(x, err, bins, errfunc=cdf_gauss, impossiblity_thresehold=1E-4):
    # [bin, x]
    probabilities = np.zeros((len(bins)-1, len(x)))
    for k in range(len(x)):
        for n in range(len(bins)-1):
            probabilities[n, k] = errfunc((bins[n + 1]-x[k])/err[k]) - errfunc((bins[n]-x[k])/err[k])
            if probabilities[n, k] < impossiblity_thresehold:
                probabilities[n, k] = 0
            elif probabilities[n, k] > (1 - impossiblity_thresehold):
                probabilities[n, k] = 1
    return probabilities


def fancy_bin(x, err, bins, errfunc=cdf_gauss, impossiblity_thresehold=1E-4):
    prob_matrix = build_prob_matrix(x, err, bins, errfunc=errfunc, impossiblity_thresehold=impossiblity_thresehold)
    hist = np.zeros((len(bins)-1, 2))
    for i in range(prob_matrix.shape[0]):
        mean, e_x_sqared = calc_mean_var_fast(prob_matrix[i])
        hist[i, 0] = mean
        hist[i, 1] = (e_x_sqared - mean**2)**0.5
    return hist


if __name__ == '__main__':
    e, e_x_sqared = calc_mean_var([0.5, 0.5, 0.5, 0.5])
    e, e_x_sqared = calc_mean_var_fast([0.5, 0.5, 0.5, 0.5])
    print(f'<x> = {e}')
    print(f'<x**2> = {e_x_sqared}')
    print(f'Var(x) = <x**2> - <x>**2 = {e_x_sqared - e**2}')

    # asinging test:
    values = np.array([4, 7, 10])
    errors = np.array([2, 4, 2])
    binning = np.linspace(0, 15, 16)

    import matplotlib.pyplot as plt
    # plt.plot(cdf_gauss(np.linspace(-4, 4, 20)))
    # plt.show()

    hist = fancy_bin(values, errors, binning)
    print(hist)
    # plt.plot(hist)
    plt.errorbar(np.linspace(0.5, 14.5, 15), hist[:, 0], yerr=hist[:, 1])
    plt.show()

    # asinging test:
    values = np.array([4, 4, 5, 6, 8, 6.5, 8, 9.5, 7, 10])
    errors = np.array([2, 4, 2, 2, 3.4, 2, 1, 2.1, 1, 1.2])
    binning = np.linspace(0, 15, 4)

    hist = fancy_bin(values, errors, binning)
    print(hist)
    plt.errorbar(np.linspace(2.5, 12.5, 3), hist[:, 0], yerr=hist[:, 1])
    plt.show()
