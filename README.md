# MTL-AMCOT-BB
MTL-AMCOT-BB

Title:Adaptive Multi-Cognitive Objective Temporal Task Approach for Predicting AD Progression   
"Adaptive Multi-Cognitive Objective Temporal Task Approach for Predicting AD Progression",has been accepted as a regular paper at IEEE International Conference on Bioinformatics and Biomedicine 2024 (IEEE BIBM 2024). 


\begin{algorithm}
    \caption{The APM-based Algorithm with Barzilai-Borwein Step Size}
    \label{alg:algorithm2}
    \begin{algorithmic}[1] % [1] enables line numbers
        \item[] \textbf{Input:}  $X, Y, \lambda_1, \lambda_2, \lambda_3$
        \item[] \textbf{Output:}  $W, \text{funcVal}$ 
        \STATE Initialization: $W_0 = W_{-1}, t_{-1} = 0, t_0 = 1, i = 1$
        \REPEAT
            \STATE $\alpha_i = \frac{t_{i-1} - 1}{t_i}$
            \STATE $S_i = W_{i-1} + \alpha_{i-1}(W_{i-1} - W_{i-2})$
            \STATE compute $s_{k-1} = W_k - W_{k-1}$, $g_{k-1} = f'(S_k) - f'(S_{k-1})$
            \STATE compute $\eta_k$ by Barzilai-Borwein Step Size
            \STATE compute $W_i = \pi(S_i - \frac{1}{\eta_k} f'(S_i))$
            \STATE $t_i = \frac{1 + \sqrt{1 + 4t_{i-1}^2}}{2}$
        \UNTIL{convergence}
    \end{algorithmic}
\end{algorithm}
