# MTL-AMCOT-BB

## Title
"Adaptive Multi-Cognitive Objective Temporal Task Approach for Predicting AD Progression" has been accepted as a regular paper at IEEE International Conference on Bioinformatics and Biomedicine 2024 (IEEE BIBM 2024).

![Target Correlation Matrix 2](/Target_correlation_5target1.png)  <!-- Assuming this is a typo and should be a different file -->
*Figure 2: Target Correlation Matrix for Target 2*

### Abstract Flowchart

Below is the flowchart describing the process used in our study:

![Flowchart Abstract](/Fig1_flowchart_abstract.png)
*Figure 3: Flowchart Abstract*

## Algorithm: APM-based Algorithm with Barzilai-Borwein Step Size

### Input:
- $X, Y, \lambda_1, \lambda_2, \lambda_3$

### Output:
- $W, \text{funcVal}$

### Procedure:
1. **Initialization:**
   - $W_0 = W_{-1}, t_{-1} = 0, t_0 = 1, i = 1$

2. **Repeat until convergence:**
   - Compute $\alpha_i = \frac{t_{i-1} - 1}{t_i}$
   - Update $S_i = W_{i-1} + \alpha_{i-1}(W_{i-1} - W_{i-2})$
   - Calculate $s_{k-1} = W_k - W_{k-1}$ and $g_{k-1} = f'(S_k) - f'(S_{k-1})$
   - Determine $\eta_k$ using the Barzilai-Borwein Step Size
   - Compute $W_i = \pi(S_i - \frac{1}{\eta_k} f'(S_i))$
   - Update $t_i = \frac{1 + \sqrt{1 + 4t_{i-1}^2}}{2}$


## Dataset ADNI
[ADNI Database](https://adni.loni.usc.edu/)

## Demographic Information of Subjects at Different Time Points

| Time Point | Attribute      | Values         |
|------------|----------------|----------------|
| **Baseline (M00)** | Sample Size     | 1532           |
|                    | (AD, MCI, CN)   | (315, 817, 400)|
|                    | (Female, Male)  | (688, 844)     |
|                    | APOE4 (0, 1, 2) | (795, 567, 160)|
| **M06**            | Sample Size     | 1317           |
|                    | (AD, MCI, CN)   | (265, 669, 383)|
|                    | (Female, Male)  | (570, 747)     |
|                    | APOE4 (0, 1, 2) | (675, 493, 149)|
| **M12**            | Sample Size     | 1375           |
|                    | (AD, MCI, CN)   | (366, 768, 241)|
|                    | (Female, Male)  | (606, 769)     |
|                    | APOE4 (0, 1, 2) | (721, 508, 146)|
| **M24**            | Sample Size     | 1099           |
|                    | (AD, MCI, CN)   | (132, 653, 314)|
|                    | (Female, Male)  | (487, 612)     |
|                    | APOE4 (0, 1, 2) | (597, 399, 103)|
| **M36**            | Sample Size     | 432            |
|                    | (AD, MCI, CN)   | (0, 277, 155)  |
|                    | (Female, Male)  | (169, 263)     |
|                    | APOE4 (0, 1, 2) | (259, 142, 31) |
| **M48**            | Sample Size     | 289            |
|                    | (AD, MCI, CN)   | (0, 172, 117)  |
|                    | (Female, Male)  | (127, 162)     |
|                    | APOE4 (0, 1, 2) | (175, 92, 22)  |
