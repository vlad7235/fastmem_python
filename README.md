# A python demonstrator of the Harmonic space-based MEM algorithm

This is a simple demonstrator of the harmonic space-based Maximum Entropy algorithm
described in detail in Hobson et al (1998) [1] and Stolyarov et al (2002) [2].

## Toy model separation
The jupyter notebook forsep_test.ipynb illustrates the method using a simple 2-component and 3-frequency model. 
All the files are provided in ./data folder. Two toy components (some arbitrary images) are stored in testcomp.fits.
The noise maps used in data construction are stored in testnoise.fits. A mixture of this two components with added noise 
is in testdata.fits. The algorithm uses only the statistical properties of the noise and components, which are assumed to be
known apriori, as well as the known mixing matrix, to find the best possible solution for component reconstruction. 
The beam in this example is assumed to be very small.

## Cosmological HI extraction
The notebook forsep_test-CosmoHI.ipynb illustrates how this approach can be used to recover non-diffuse signal poorly correlated
between the frequency channels. The full data model includes 2048 frequency maps with synchrotron, free-free and cosmological HI
contribution at the frequencies from 88 to 202 MHz for some particular EoR model taken from 21SSD database [3]. The synchrotron
is assumed to have a varying spectral index.

To extract the diffuse component(s) the full dataset is not requred. Instead, the preconditioning of the data is performed, making one frequency
map from the 64 neighbouring channels. This reduced dataset, containing only 32 frequency maps, is used to extract the diffuse components - one mean synchrotron,
and another one accomodates the synchrotron variations due to the variable spectral index, and free-free.

After the diffuse components are recovered, we can model their contribution at each of 2048 frequencies, and subtract it from the initial dataset.
This would give us the difference maps containing cosmological HI, instrumental noise and the residuals from the diffuse components.

In this particular example the instrumental noise is assumed to be very small (long integration limit), the beam is also small (SKA-MID case) and the statistical
properties of the cosmological HI are assumed to be known.

Since the data files are very big, they are not uploaded to the GitHub.


## Acknowledgement

This research was supported by H2020-Astronomy ESFRI and Research
Infrastructure Cluster (Grant Agreement number: 653477).

## References

1. Hobson et al, Foreground separation methods for satellite observations of the cosmic microwave background, MNRAS, 300, p1 (1998)
2. Stolyarov et al, All-sky component separation for the Planck mission, MNRAS, 336, p97 (2002)
3. Semelin et al, 21SSD: a public data base of simulated 21-cm signals from the epoch of reionization, MNRAS, 472, p4508 (2017)
