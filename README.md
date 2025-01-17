# IsoFracPy

This project is developed by Kanon Kino (kanon@hydra.t.u-tokyo.ac.jp).
Code refactoring was mostly done by ChatGPT-4o.

## Project Structure

```
IsoFracPy/
│
├── IsotopeFractionationModel/
│   ├── __init__.py
│   ├── config.py
│   ├── BasicUtility.py
│   ├── EquilibriumFractionation.py
│   ├── KineticFractionation.py
│   ├── SeaEvaporationIsotopeCalculation.py
│   ├── InitialCondition.py
│   ├── RayleighDistillation.py
│   ├── PostPrecipitationProcess.py
│
├── main.py 
└── tests/  
```

### Dependencies

```
config.py
    ├── BasicUtility.py
    │　　　└── config.py
    ├── EquilibriumFractionation.py
    │     ├── BasicUtility.py
    │　　　└── config.py
    ├── KineticFractionation.py
    │     ├── BasicUtility.py
    │     ├── EquilibriumFractionation.py
    │　　　└── config.py
    ├── SeaEvaporationIsotopeCalculation.py
    │     ├── BasicUtility.py
    │　　　└── config.py
    ├── InitialCondition.py
    │     ├── BasicUtility.py
    │     ├── EquilibriumFractionation.py
    │     ├── KineticFractionation.py
    │     ├── SeaEvaporationIsotopeCalculation.py
    │　　　└── config.py
    ├── RayleighDistillation.py
    │     ├── BasicUtility.py
    │     ├── EquilibriumFractionation.py
    │     ├── KineticFractionation.py
    │     ├── SeaEvaporationIsotopeCalculation.py
    │     ├── InitialCondition.py
    │　　　└── config.py
    └── PostPrecipitationProcess.py
     　　　└── config.py
```