![Visits Badge](https://badges.pufler.dev/visits/WUST-FOG/cgnlse-python)[![DOI](https://zenodo.org/badge/477591470.svg)](https://zenodo.org/badge/latestdoi/477591470)

# Nonlinear phenomena in birefringent microstructured fibers

We implemented model based on two coupled nonlinear Schrödinger equations that include both the Raman and the Kerr nonlinearities. We used it to study evolution of nonlinear phenomena in the temporal and spectral domains in optical fibers exhibiting high and low birefringence.

![soliton_traping](./data/191119_polarisation_21m_lambda_1560nm_power_41mW.png)

## Usage

### Installation

1. Create a virtual environment with `python -m venv cgnlse` or using `conda` by `conda create -n cgnlse python=3.8`.
2. Activate it with `source cgnlse/bin/activate` or `conda activate cgnlse`.
3. Install `gnlse` package `pip install gnlse==2.0.0`
3. Clone this repository `git clone https://github.com/WUST-FOG/cgnlse-python.git`

```bash
python -m venv cgnlse
source cnlse/bin/activate
pip install -r requirements.txt
git clone https://github.com/WUST-FOG/cgnlse-python.git
cd cgnlse-python
```

### Examples

**Soliton trapping and orthogonal Raman scattering**

Run test script to generate above figure and reproduce the manuscript results:

```bash
python draw_soliton_traping.py
```

Note that we also provided script tu run simulations (`run_soliton_traping.py`),
however used input data is not publicly available at this time,
but may be obtained from the authors upon reasonable request.

_Inspiration: K. Stefańska et al., Soliton trapping and orthogonal Raman scattering in a birefringent microstructured fiber_

**Modulation instability in highly birefringent fibers**

To run example of vector modulation instability in highly birefringent
fibers with circularly polarized modes in the normal dispersion regime
type:

```bash
python run_modulation_instability.py
```

Note that using also raman_polarisation and setting solver.fr to 0
one can simulate the case of low-birefringent fibers.

_Inspiration: K. Zołnacz et al., Vector modulation instability in highly birefringent fibers with circularly polarized eigenmodes_

## Acknowledgement

cnlse-python is a Python set of scripts for solving
Coupled Nonlinear Schrodringer Equation. It is one of the WUST-FOG
projects developed by [Fiber Optics Group, WUST](http://www.fog.pwr.edu.pl/).

The python code based on `gnlse` package, available at 
[https://github.com/WUST-FOG/gnlse-python](https://github.com/WUST-FOG/gnlse-python).

## Citation
If you find this code useful in your research, please consider citing:

[Soliton trapping and orthogonal Raman scattering in a birefringent photonic crystal fiber](https://arxiv.org/abs/2204.13773v1):

```
@article{Stefanska:22,
      author = {Karolina Stefa\'{n}ska and Sylwia Majchrowska and Karolina Gemza
                and Grzegorz Sobo\'{n} and Jaros{\l}aw Sotor and Pawe{\l} Mergo
                and Karol Tarnowski and Tadeusz Martynkien},
      journal = {Opt. Lett.},
      number = {16},
      pages = {4183--4186},
      publisher = {Optica Publishing Group},
      title = {Soliton trapping and orthogonal Raman scattering
               in a birefringent photonic crystal fiber},
      volume = {47},
      month = {Aug},
      year = {2022},
      url = {http://opg.optica.org/ol/abstract.cfm?URI=ol-47-16-4183},
      doi = {10.1364/OL.463643}
}
```

[gnlse-python: Open Source Software to Simulate Nonlinear Light Propagation In Optical Fibers](https://arxiv.org/abs/2110.00298):

```
@misc{redman2021gnlsepython,
      title={gnlse-python: Open Source Software to Simulate
             Nonlinear Light Propagation In Optical Fibers}, 
      author={Pawel Redman and Magdalena Zatorska and Adam Pawlowski
              and Daniel Szulc and Sylwia Majchrowska and Karol Tarnowski},
      year={2021},
      eprint={2110.00298},
      archivePrefix={arXiv},
      primaryClass={physics.optics}
}
```

## License
[MIT](https://choosealicense.com/licenses/mit/)
