Scores Command
==============

The ``scores`` command calculates fragmentomics scores like the Windowed Protection Score (WPS) and motif diversity metrics.

Syntax
------

.. code-block:: bash

   pyfrag.py scores -o <OUT_DIR> [OPTIONS]

Options
-------

.. code-block:: bash

   --frag-file PATH       Single fragment file to analyze
   --frag-dir PATH        Directory containing fragment files
   --bed-file PATH        BED file defining genomic regions (required)

Examples
--------

.. code-block:: bash

   pyfrag.py --out-dir scores/ scores --frag-file sample.frag --bed-file regions.bed


BED File Format
---------------

The BED file should contain genomic regions of interest for which the windowed protection score is then calculated:

.. code-block:: text

   chr1    1000000    1001000    region1
   chr1    2000000    2001000    region2
   chr2    500000     501000     region3

Required columns:
- Column 1: Chromosome name
- Column 2: Start position (0-based)
- Column 3: End position (0-based)
- Column 4: Region name (optional)

Output
------

As output, this subcommand saves the following *CSV* files to disk:

* ``wps_<sample>.csv`` — per-position Windowed Protection Score for all regions in the input *BED* file. The file also has a ``wps_smooth`` column containing the coverage-normalised, detrended, Savitzky-Golay-smoothed signal that is suitable for between-sample comparison.
* ``wps_metrics_<sample>.csv`` — per-region WPS summary metrics (number of peaks, median peak spacing, median peak prominence and width, and the 170 bp spectral power ``power_170``). These are also used by the ``differential_wps.py`` utility script.
* ``global_sample_metrics.csv`` — one row per sample, providing end-motif diversity (Shannon entropy and Simpson index for 5' and 3' motifs) together with cohort-level WPS summaries (median ``power_170``, median peaks per region, median peak spacing, median peak prominence).

A genome-wide WPS line plot is additionally written per sample as a visual summary, together with a diagnostic ``<sample>_wps_spectrum.png`` that shows the per-region magnitude spectra and the mean curve across regions with the nucleosome-repeat band highlighted.

Calculated Metrics
------------------

Windowed Protection Score (WPS)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For every genomic position :math:`i`, the WPS measures the depletion of fragment ends within a window surrounding that position:

.. math::

   \text{WPS}(i) = N_{span}(i) - N_{end}(i)



Large values indicate protection (fewer fragment ends), whereas smaller values indicate chromatin accessibility (more fragment ends). Because this metric is not independent of local fragment depth, the output *CSV* file contains positional fragment depth, the number of ending, and the number of spanning fragments as additional fields. This way, users can calculate fractions instead of differences if that makes more sense for their datasets.

Smoothed WPS and per-region metrics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As pointed out above, raw WPS amplitudes are dominated by sequencing depth and by long-wavelength coverage drift, which make between-sample comparisons difficult. To robustly isolate nucleosome phasing signals, the ``scores`` subcommand runs the smoothing pipeline implemented in :func:`pyfraglib.scores.smooth_wps`:

1. coverage-normalise via ``wps / (median_depth + 1)`` per region,
2. subtract a 1001 bp rolling median,
3. apply a Savitzky-Golay filter (window 21 bp, polynomial order 3).

The resulting ``wps_smooth`` column is appended to ``wps_<sample>.csv``. :func:`pyfraglib.scores.wps_region_metrics` then derives a set of interpretable per-region metrics from the same pipeline and writes them to ``wps_metrics_<sample>.csv``:

* ``n_peaks`` — number of peaks detected with ``scipy.signal.find_peaks`` using a prominence threshold of 0.02 x per-region IQR and a minimum peak-to-peak distance of 120 bp.
* ``mean_peak_dist`` / ``median_peak_dist`` — mean and median of consecutive peak-to-peak spacings within the region.
* ``median_prom`` / ``median_width`` — median peak prominence and median peak width.
* ``power_170`` — relative spectral power in the nucleosome-repeat band (see below).
* ``n_positions`` — number of input positions for the region.

Parameter choices
~~~~~~~~~~~~~~~~~

The three smoothing parameters are tied directly to the average nucleosome-repeat period :math:`\lambda = 170` bp, which sets the target frequency :math:`f_\lambda = 1/\lambda \approx 0.0059` cycles/bp. Each parameter is chosen so that :math:`f_\lambda` sits comfortably inside the passband of the combined filter.

*1001 bp rolling-median detrend window.* Rolling-median subtraction acts as a high-pass filter with effective cutoff near :math:`1/W` cycles/bp. With :math:`W = 1001`, the cutoff lies roughly a factor of six below :math:`f_\lambda`, so the nucleosome oscillation passes essentially unattenuated while long-wavelength coverage drift (kilobase-scale protected / unprotected blocks that otherwise dominate raw WPS amplitude) is removed. The window is adaptively shrunk for short regions via :math:`\min(W,\; \max(101,\; (L/3)\,|\,1))`.

*21 bp Savitzky-Golay window.* A Savitzky-Golay filter with window :math:`N` has its 3 dB point near :math:`1/(N/2)` cycles/bp; at :math:`N = 21` that corresponds to a period floor of roughly 10 bp -- more than an order of magnitude below :math:`\lambda`, so the 170 bp oscillation is preserved while per-base Poisson noise in the spanning-minus-ending count is averaged out.

Spectral analysis of WPS
~~~~~~~~~~~~~~~~~~~~~~~~

Let :math:`x[n]` denote the detrended, coverage-normalised WPS at positions :math:`n = 0, 1, \ldots, N - 1` where :math:`N` is the size of the region. We first remove the zero-frequency component,

.. math::

   \tilde{x}[n] = x[n] - \frac{1}{N} \sum_{k=0}^{N-1} x[k],

then compute the one-sided discrete Fourier transform on a zero-padded signal of length :math:`N_{\text{fft}} = 2^{\lceil \log_2 N \rceil}` using the ``numpy.fft.rfft`` implementation:

.. math::

   X[k] = \sum_{n=0}^{N_{\text{fft}} - 1} \tilde{x}[n]\, e^{-j 2\pi k n / N_{\text{fft}}},
   \qquad k = 0, 1, \ldots, N_{\text{fft}} / 2.

The corresponding frequency grid is :math:`f_k = k / N_{\text{fft}}`, expressed in cycles per base pair because the per-position sampling interval is 1 bp. The magnitude spectrum :math:`S[k] = |X[k]|` is then integrated over a narrow band around the nucleosome repeat length :math:`\lambda = 170` bp. With

.. math::

   \mathcal{B} = \left\{\, k : \frac{1}{\lambda + 10} < f_k < \frac{1}{\lambda - 10} \,\right\},

the 170 bp spectral-power feature is the in-band mean magnitude normalised by the overall mean magnitude,

.. math::

   \texttt{power\_170}
       = \frac{\displaystyle \frac{1}{|\mathcal{B}|} \sum_{k \in \mathcal{B}} S[k]}
              {\displaystyle \frac{1}{N_{\text{fft}}/2 + 1} \sum_{k=0}^{N_{\text{fft}}/2} S[k]}.


This ratio is dimensionless and independent of depth-driven amplitude. The mean-over-bins normalisation is chosen so that a flat spectrum (i.e. no preferred periodicity) yields numerator and denominator of equal magnitude, giving power_170 = 1. Values above 1 indicate that the nucleosome-repeat band carries more energy per bin than the spectrum-wide average (i.e. nucleosomes are regularly phased).

Motif Diversity
~~~~~~~~~~~~~~~

Given an end motif distribution :math:`X = (x_1, ..., x_n)` where :math:`x_i` is the observed count for end motif :math:`i`, we calculate Shannon entropy :math:`H(X)` as

.. math::

   H(X) = -\sum_{i=1}^{n} p(x_i) \log_2 p(x_i)

and Simpson index :math:`D(X)` as

.. math::

    D(X) = \sum_{i=1}^{n} x_i^2

.
