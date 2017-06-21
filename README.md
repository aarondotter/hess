# Python for Hess Diagrams

## Requirements

hess.py requires:
* [read_mist_models](https://github.com/jieunchoi/MIST_codes) for isochrone manipulation
* [shapely](https://github.com/Toblerity/Shapely) for geometry.

## Usage

Read in MIST isochrone. The code is set up to use the Bessell V,V-I CMD by default.
```python
p000=ISOCMD('MIST_v1.0_feh_p0.00_afe_p0.0_vvcrit0.4_UBVRIplus.iso.cmd')

#isochrone logAge=9.6, excluding TP- and post-AGB
iso=p000.isocmds[92]
iso=iso[iso['EEP']<=808]
```
Create a rectangular Hess Diagram, attach an isochrone, and calculate IMF weights.

```python
H=HessDiagram(0.5,8,2.5,23,100,100)
H.add_iso(iso,DM=10)
H.do_mass_intervals()
H.do_IMF(alpha=2.35)
```

