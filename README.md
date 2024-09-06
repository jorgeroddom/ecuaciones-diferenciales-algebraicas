En este repositorio vamos a aproximar numéricamente las soluciones de sistemas de ecuaciones diferenciales algebraicas de índice 1, es decir, problemas de la forma:
```math
\begin{cases}
    y'= f(y,z), \\
    0 = g(y,z),
\end{cases}
```
donde $f: \mathbb{R}^n \times \mathbb{R}^m \rightarrow \mathbb{R}^n$, $g: \mathbb{R}^n \times \mathbb{R}^m \rightarrow \mathbb{R}^m$ y y $g_z$ es invertible en un entorno de la solución, siendo
```math
g_z(y,z) = \begin{bmatrix}
\dfrac{\partial g_1 }{\partial z_1}(y,z) & \ldots & \dfrac{\partial g_1 }{\partial z_m}(y,z)\\
\vdots & \ddots & \vdots\\
\dfrac{\partial g_m }{\partial z_1}(y,z) & \ldots & \dfrac{\partial g_m }{\partial z_m}(y,z)
\end{bmatrix}.
```
Dependiendo del valor de $n$ y $m$ distinguimos los siguientes archivos:
- $n=1$ y $m=1$: `state space form n1m1.py` y `epsilon encajado n1m1.py`.
- $n \ge 2$ y $m=1$: `state space form n2m1.py`.
- $n \ge 2$ y $m \ge 2$: `state space form n2m2.py` y `epsilon encajado n2m2.py`.

También se pueden consultar las funciones y tableros de Butcher utilizados durante mi TFG:
- Tableros de Butcher: `tableros de butcher.py`.
- Funciones péndulo simple: `pendulo simple.py`.
- Funciones amplificador de transistores: `amplificador de transistores.py`.

Para más información, se puede consultar `Resumen.pdf`.
