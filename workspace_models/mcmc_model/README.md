## Setup
### OS
We have verified this code on only MacOSX(10.14.6)
### The version of Python
- Only Python3
- 3.6 or upper
### Required packages
- PyStan
- numpy
- matplotlib
- seaborn
- pandas
- optuna

### Run optimization code
Before all, create directories `./data/csv`, `./data/img`, `./data/json`, and `psycho_ex`.
#### MCMC model
Run `python optimize_mcmc.py`

#### EKF model
Run `python optimize_ekf.py`

### Run test code
#### MCMC model
Run `python test_mcmc.py`

#### EKF model
Run `python test_ekf.py`

#### Visualizing the results
You can visualize the results with The [Jupyter notebook](https://jupyter.org/). Run `jupyter notebook` or `jupyter lab`, and open `visualize.ipynb`.

## References
### Extended Kalman Filter
- [Wikipedia page](https://en.wikipedia.org/wiki/Extended_Kalman_filter)
- [Post on Towards Data Science](https://towardsdatascience.com/extended-kalman-filter-43e52b16757d)
- [Cognitive/Neuroscientific perspective of Kalman Filter](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4415408/)
### Change Blindness
- [Model of prediction-related change blindness](https://books.google.co.uk/books?id=clVkDQAAQBAJ&pg=PA387&lpg=PA387&dq=NEURAL+MODELS+OF+PREDICTION+AND+SUSTAINED+INATTENTIONAL+BLINDNESS&source=bl&ots=9x_pkVTxli&sig=ACfU3U2uQbIYlzfovTGK-jdmACCKcbrdqg&hl=en&sa=X&ved=2ahUKEwi3g9-ojcXgAhXBXhUIHQpJAGMQ6AEwAnoECAcQAQ#v=onepage&q=NEURAL%20MODELS%20OF%20PREDICTION%20AND%20SUSTAINED%20INATTENTIONAL%20BLINDNESS&f=false)


## Models
### MCMC
We use [Stan](https://mc-stan.org/), which is a platform for statistical modeling and high-performance statistical computation, especially approximative Bayesian estimation. The overview of probablistic models is written on `.stan` files. The inference is done by the following procedures.
- Import PyStan `import pystan`
- Define the model `sm = pystan.StanModel(file='path/to/StanFile')`
- Run the inference (MCMC) `sm.sampling(data=stan_data, iter=5000, chains=1)`

### EKF (Extended Kalman Filter)
This involves two phases (See `optimize_ekf.py`)
- `predict_phase`
- `update_phase`