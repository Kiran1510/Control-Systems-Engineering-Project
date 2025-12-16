# Control Systems Engineering Project
**Anti-Lock Braking System (ABS) Controller Design and Analysis**

This repository contains comprehensive analysis and implementation of observer-based state feedback control for automotive anti-lock braking systems using quarter-vehicle dynamics models.

## Author
**Kiran Sairam Bethi Balagangadaran**  
MS Robotics, Northeastern University  
Professor: Dr. Nivii Kalavakonda Chandrasekar  
Course: ME5659 - Control Systems Engineering

## Project Overview

This project demonstrates advanced control systems design techniques applied to automotive safety systems. Two comprehensive reports are included:

### 1. [Final Report - Observer-Based ABS Controller](Final%20Report%20Controls%20Systems.pdf)
**Complete system analysis and controller design for ABS**
- Nonlinear quarter-vehicle model with tire-road friction dynamics
- Linearization around optimal slip ratio (λ* = 0.2, Vx* = 20 m/s)
- Stability, controllability, and observability analysis
- Standard form decomposition revealing unobservable position state
- State feedback controller design using pole placement
- Luenberger observer design for state estimation from slip measurements
- Validation on nonlinear plant with disturbance rejection
- **Key Results:** 50% overshoot, 0.039s settling time, excellent disturbance rejection

### 2. [Extra Credit - Advanced Control Techniques](Extra%20Credit%20Control%20Systems.pdf)
**Minimum energy control and animated visualization**
- Minimum energy open-loop steering to equilibrium
- Observer convergence demonstration with poor initial estimates
- Real-time animated simulation showing ABS operation
- Interactive visualization with heads-up display (HUD)
- **Key Results:** 0.5s observer convergence, successful regulation with animation

## Technical Highlights

### System Model
- **Plant:** Nonlinear quarter-vehicle ABS model
- **States:** Stopping distance (Sx), vehicle velocity (Vx), slip ratio (λ)
- **Control Input:** Brake torque (0-1200 Nm)
- **Output:** Slip ratio from wheel speed sensors
- **Operating Point:** λ* = 0.2, Vx* = 20 m/s (maximum friction coefficient)

### Control Architecture
```
Wheel Speed Sensors → Observer → State Feedback Controller → Brake Modulator → Vehicle Dynamics
                        ↑                                           ↓
                        └─────────── Slip Ratio Measurement ────────┘
```

### Key Achievements
- **Stability Analysis:** Identified unstable open-loop system with eigenvalues at {0, 0.0203, 2.9059}
- **Controllability:** Full rank controllability matrix - all states controllable via brake torque
- **Observability:** Rank-2 observable subsystem [Vx, λ] sufficient for ABS control
- **Controller Design:** Pole placement with ts = 2s, ζ = 0.9 specifications
- **Observer Design:** 5× faster poles than controller for rapid convergence (0.446s)
- **Performance:** 50% overshoot, 0.039s settling time on nonlinear plant
- **Disturbance Rejection:** Maintains performance under 45% friction reduction (dry→wet)
- **Minimum Energy Control:** Optimal open-loop steering within actuator limits
- **Animation:** Real-time visualization demonstrating control effectiveness

## Repository Structure
```
Control-Systems-Engineering-Project/
├── Matlab Code/           # MATLAB simulation files
│   ├── sec41.m           # Section 4.1: State feedback controller
│   ├── sec42.m           # Section 4.2: Luenberger observer
│   ├── sec43.m           # Section 4.3: Observer-based control (linear)
│   ├── sec44.m           # Section 4.4: Nonlinear plant application
│   ├── sec45.m           # Section 4.5: Disturbance rejection
│   ├── sec51.m           # Section 5.1: Minimum energy control
│   ├── sec52.m           # Section 5.2: Observer convergence demo
│   └── sec53.m           # Section 5.3: Animated ABS simulation
│   └── sec53.gif         # Animation output (6.74 MB)
├── Extra Credit Control Systems.pdf       # Extra credit report
├── Final Report Controls Systems.pdf      # Main project report
└── README.md                              # This file
```

## Visualizations

### System Animation
The project includes an animated simulation (`sec53.gif`) showing:
- Vehicle position along roadway
- Real-time velocity tracking (true vs. estimated)
- Real-time slip ratio tracking (true vs. estimated)
- Heads-up display with system metrics
- Observer convergence from poor initial estimates

### Static Plots
All sections include comprehensive plots:
- State trajectories with equilibrium tracking
- Control effort within actuator limits
- Phase portraits showing system dynamics
- Estimation error decay (exponential convergence)
- Disturbance rejection response

## Getting Started

### Prerequisites
- MATLAB R2020b or later
- Control System Toolbox
- (Optional) Symbolic Math Toolbox for analytical verification

### Running the Code

Each MATLAB file is self-contained and can be run independently:
```matlab
% State feedback controller
run('sec41.m')

% Luenberger observer
run('sec42.m')

% Observer-based control on linear plant
run('sec43.m')

% Nonlinear plant validation
run('sec44.m')

% Disturbance rejection
run('sec45.m')

% Minimum energy control
run('sec51.m')

% Observer convergence demonstration
run('sec52.m')

% Animated ABS simulation
run('sec53.m')
```

### Animation Instructions

To view the ABS animation:
1. Open `sec53.gif` in any image viewer
2. Or run `sec53.m` in MATLAB to generate a new animation
3. Online version available at: [Google Drive Link](https://drive.google.com/file/d/1JR51Bv2CDtsQCkMc6gjYQltOd6ImC8ex/view?usp=sharing)

## Technical Details

### Mathematical Formulation

**Nonlinear State Equations:**
```
ẋ₁ = x₂
ẋ₂ = -μ(x₃,x₂)Fₙ/m
ẋ₃ = -μ(x₃,x₂)Fₙ/x₂[(1-x₃)/m + R²/Jw] + R/(Jwx₂)u
```

**Friction Model (Pacejka-style):**
```
μ(λ,Vx) = c₁(1-e^(-c₂λ)) - c₃λ · e^(-c₄Vx)
```

**Linearized System:**
```
A = [0      1.0000   0     ]    B = [0      ]
    [0      0.1883   1.4362]        [0      ]
    [0      0.3178   2.7380]        [0.0146 ]

C = [0  0  1]
```

**Controller Design:**
- Desired poles: p = -2.000 ± 0.969j
- Feedback gains: K = [294.887, 474.404]

**Observer Design:**
- Observer poles: pobs = -10.000 ± 4.843j (5× faster)
- Observer gains: L = [401.871; 22.926]

### Performance Metrics

| Metric | Linear Plant | Nonlinear Plant |
|--------|-------------|-----------------|
| Overshoot | 171.36% | 50.17% |
| Settling Time | 2.907s | 0.039s |
| Steady-State Error | 0.000128 | 0.000000 |
| Control Effort | 656-765 Nm | 0-725 Nm |

### Vehicle Parameters

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Vehicle Mass | m | 342 kg |
| Wheel Inertia | Jw | 1.13 kg·m² |
| Wheel Radius | R | 0.33 m |
| Normal Force | Fₙ | 3354 N |
| Gravity | g | 9.81 m/s² |

## Key Findings

### 1. Stability Analysis
- Open-loop system is **unstable** with two positive eigenvalues
- Instability arises from operating on downslope of friction curve (∂μ/∂λ < 0)
- Active feedback control is **necessary** for ABS operation

### 2. Controllability/Observability
- System is **fully controllable** via brake torque
- Position state is **unobservable** from slip measurements
- Observable subsystem [Vx, λ] is **sufficient** for ABS objectives

### 3. Controller Performance
- **Separation principle validated:** Observer-based control nearly identical to perfect state feedback
- **Nonlinear improvement:** 121% overshoot reduction due to friction saturation
- **Disturbance rejection:** Maintains performance under 45% friction drop

### 4. Practical Implications
- Controller uses **only wheel speed sensors** (no direct velocity measurement needed)
- **Bounded control effort** within actuator limits (0-1200 Nm)
- **Fast convergence:** Observer estimates ready within 0.5s
- **Robust to model mismatch:** Performs better on nonlinear plant than predicted

## Applications

This control strategy is applicable to:
- **Automotive Safety Systems:** Modern ABS implementations
- **Brake-by-Wire Systems:** Electronic brake force distribution
- **Autonomous Vehicles:** Emergency braking controllers
- **Motorcycle ABS:** Similar slip ratio control objectives
- **Aircraft Landing Systems:** Anti-skid braking control

## Future Work

### Proposed Extensions
1. **LQR Optimization:** Replace pole placement with Linear Quadratic Regulator for optimal trade-off between performance and control effort
2. **Adaptive Observer:** Extended Kalman Filter to handle model mismatch on nonlinear plant
3. **Full Vehicle Model:** Extend to 4-wheel coordination with load transfer dynamics
4. **Road Surface Estimation:** Adaptive slip reference based on detected friction coefficient
5. **Multi-Objective Control:** Simultaneously optimize stopping distance and ride comfort

### Unmodeled Effects to Include
- Load transfer (front/rear weight distribution during braking)
- Road grade compensation (uphill/downhill operation)
- Tire temperature and wear dynamics
- Aerodynamic drag forces
- Suspension dynamics coupling

## References

### Control Theory
1. K. Ogata, *Modern Control Engineering*, 5th ed. Prentice Hall, 2010
2. H. K. Khalil, *Nonlinear Systems*, 3rd ed. Prentice Hall, 2002
3. B. Friedland, *Control System Design: An Introduction to State-Space Methods*, Dover, 2005

### Automotive Engineering
4. U. Kiencke, L. Nielsen, *Automotive Control Systems*, 2nd ed. Springer, 2004
5. R. Isermann, *Automotive Control*, 1st ed. Springer, 2021
6. H. Pacejka, *Tire and Vehicle Dynamics*, 3rd ed. Elsevier, 2012
7. R. Rajamani, *Vehicle Dynamics and Control*, 2nd ed. Springer, 2012

### ABS-Specific Literature
8. P. B. Bhivate, "Modelling & development of antilock braking system," B.Tech Thesis, NIT Rourkela, 2011
9. A. B. Sharkawy, "Genetic fuzzy self-tuning PID controllers for ABS," *Eng. App. of AI*, vol. 23, pp. 1041-1052, 2010
10. H. Mirzaeinejad, M. Mirzaei, "Non-linear control of wheel slip in ABS," *Control Eng. Practice*, vol. 18, pp. 918-926, 2010
11. S. Ç. Başlamişli et al., "Robust control of ABS," *Vehicle System Dynamics*, vol. 45, no. 3, pp. 217-232, 2007
12. S. B. Choi, "ABS with continuous wheel slip control," *IEEE Trans. Control Syst. Tech.*, vol. 16, no. 5, 2008

## Tools Used

- **MATLAB/Simulink** - All analysis, simulations, and computations
- **LaTeX/Overleaf** - Professional technical report preparation
- **Claude AI** - Document formatting and structuring assistance
- **Microsoft Copilot** - MATLAB debugging and code optimization

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

Special thanks to **Professor Dr. Nivii Kalavakonda Chandrasekar** for guidance throughout this project and for her expertise in control systems engineering at Northeastern University.

## Contact

**Kiran Sairam Bethi Balagangadaran**  
MS Robotics, Northeastern University  
Email: bethi.k@northeastern.edu  
GitHub: [@Kiran1510](https://github.com/Kiran1510)

---

*This project demonstrates advanced control system design techniques including linearization, stability analysis, controllability/observability decomposition, state feedback control, Luenberger observer design, and validation on nonlinear dynamics - all applied to a practical automotive safety system.*
