```mermaid
graph TD
    %% Data sources and Python processing layer
    subgraph Layer_1 [1. Pre-processing Pipeline]
        Raw[(data/raw/)] -->|python/ scripts| Deriv[(data/derived/)]
    end

    %% Main R analysis workflow (3 Chapters)
    subgraph Layer_2 [2. Analysis Pipeline]
        Deriv -->|scripts/chap1_exploration/| Env{R Environment}
        Env -->|scripts/chap2_modelling/| Env
        Env -->|scripts/chap3_predictions/| Out([outputs/])
    end

    %% R functions, testing, and documentation
    subgraph Layer_3 [3. Tools, Validation & Docs]
        Fcts[R/ functions] -.->|functions_loader.R| Env
        Tests[tests/testthat/] -.->|Validates| Fcts
        Vignettes[vignettes/] -.->|Explains| Fcts
    end
