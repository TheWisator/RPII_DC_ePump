import optuna
import numpy as np

from pump_design_prop import *


def objective(trial):
    n = trial.suggest_int('n', 35_000, 38_001)
    psi = trial.suggest_float('psi', 0.69, 0.705)
    phi = trial.suggest_float('phi', 0.145, 0.155)
    R = trial.suggest_float('R', 0.6640, 0.6650)
    delta_h = 0.25
    
    
    hyperparameters = {
        'n': n,
        'psi':psi,
        'phi':phi,
        'R':R,
        'delta_h':delta_h,
    }

    accuracy = run(hyperparameters)

    # Handle pruning based on the intermediate value.    
    if trial.should_prune():
        raise optuna.exceptions.TrialPruned()
    
    return accuracy


if __name__ == "__main__":
    study = optuna.create_study(
        direction='maximize',
        storage="sqlite:///db.sqlite3"
        )
    study.optimize(objective, n_trials=1_000_000)  # You can adjust the number of trials


    pruned_trials = study.get_trials(states=(optuna.trial.TrialState.PRUNED,))
    complete_trials = study.get_trials(states=(optuna.trial.TrialState.COMPLETE,))

    print("Study statistics: ")
    print("  Number of finished trials: ", len(study.trials))
    print("  Number of pruned trials: ", len(pruned_trials))
    print("  Number of complete trials: ", len(complete_trials))

    # Print the best hyperparameters
    #print('Best hyperparameters:', study.best_params)

    # Retrieve the best hyperparameters
    #best_hyperparams = study.best_params

    # Create an instance of your TD3 model with the best hyperparameters
    #run(best_hyperparams, final=True)
