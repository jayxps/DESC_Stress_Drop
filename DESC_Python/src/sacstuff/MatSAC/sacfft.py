import numpy as np


def sacfft(data, tt, ifplot=False):
    """
    Translate time-domain data to frequency domain.
    Returns (MX, Phase, ff) analogous to sacfft.m.
    """
    if len(tt) < 2:
        raise ValueError("Time vector must have at least 2 points")
    delta = float(tt[1] - tt[0])
    Fs = 1.0 / delta
    Fn = Fs / 2.0

    N = len(data)
    NFFT = 1 << int(np.ceil(np.log2(N)))

    FFTX = np.fft.fft(data, NFFT)
    NumUniquePts = int(np.ceil((NFFT + 1) / 2.0))
    FFTX = FFTX[:NumUniquePts]

    MX = np.abs(FFTX)
    Phase = np.angle(FFTX)

    # Glenn's approach scaling
    MX = MX * delta

    ff = np.arange(NumUniquePts) * 2.0 * Fn / NFFT

    if ifplot:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.subplot(2, 1, 1)
        plt.plot(tt, data)
        plt.xlabel('Time')
        plt.title('Data in the Time Domain')
        plt.subplot(2, 1, 2)
        plt.loglog(ff, MX)
        plt.xlabel('Frequency')
        plt.title('Data in the Frequency Domain')
        plt.tight_layout()
        plt.show()

    return MX, Phase, ff
