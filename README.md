# cnlse-python

cnlse-python is a Python set of scripts for solving
Coupled Nonlinear Schrodringer Equation. It is one of the WUST-FOG students
projects developed by [Fiber Optics Group, WUST](http://www.fog.pwr.edu.pl/).

![soliton_traping](./data/191119_polarisation_21m_lambda_1560nm_power_41mW.png)

## Installation

1. Create a virtual environment with `python -m venv cnlse` or using `conda` by `conda create -n cnlse python=3.8`.
2. Activate it with `. cnlse/bin/activate` or `conda activate cnlse`.
3. Install `gnlse-python` package `pip install git+https://github.com/WUST-FOG/gnlse-python.git`
3. Clone this repository `git clone https://github.com/WUST-FOG/vnlse-python.git`

```bash
python -m venv cnlse
. cnlse/bin/activate
pip install git+https://github.com/WUST-FOG/gnlse-python.git
git clone https://github.com/WUST-FOG/cnlse-python.git
cd gnlse-python
```

Run test script:

```bash
python run_soliton_traping.py
```

## License
[MIT](https://choosealicense.com/licenses/mit/)
