# Definitions
* Passive rotation: transforms a vector expressed in one frame to a vector expressed in another frame. Remapping of coordinates, frames are static, vector is unchanged.
* Active rotation: rotates a vector expressed in one coordinate frame; does not change the coordinate frame the vector is expressed in.

# Convention
#### Right-handed (Hamilton) product:
* $\mathbf{i}^2 = \mathbf{j}^2 = \mathbf{k}^2 = \mathbf{i}\mathbf{j}\mathbf{k} = -1$
  * $\mathbf{i}\mathbf{j} = -\mathbf{j}\mathbf{i} = \mathbf{k}$
  * $\mathbf{j}\mathbf{k} = -\mathbf{k}\mathbf{j} = \mathbf{i}$
  * $\mathbf{k}\mathbf{i} = -\mathbf{i}\mathbf{k} = \mathbf{j}$
 
#### Scalar-first vector form:
  * $ \mathbf{q} = \begin{bmatrix} q_s \\ \mathbf{q}_v \end{bmatrix} = \begin{bmatrix} q_s \\ q_x \\ q_y \\ q_z \end{bmatrix} $

#### Quaternion product: $\odot$
* $\mathbf{p} \odot \mathbf{q} = L(\mathbf{p}) \mathbf{q} = R(\mathbf{q})\mathbf{p}$
  * $
L(\mathbf{p}) = 
\begin{bmatrix}
p_s & -\mathbf{p}_v^T \\
\mathbf{p}_v & p_s I +[\mathbf{p}_v]_\times
\end{bmatrix}
$
  * $
R(\mathbf{q}) = 
\begin{bmatrix}
q_s & -\mathbf{q}_v^T \\
\mathbf{q}_v & q_s I -[\mathbf{q}_v]_\times
\end{bmatrix}
$
  * Because left- and right-multiplication commute, $L(\mathbf{p})$ and $R(\mathbf{q})$ always commute:
$
R(\mathbf{r})L(\mathbf{p}) \mathbf{q} = (\mathbf{p} ⊙ \mathbf{q}) ⊙ \mathbf{r}) = \mathbf{p} ⊙ (\mathbf{q} ⊙ \mathbf{r}) = L(\mathbf{p}) R(\mathbf{r}) \mathbf{q}
$
 
#### Conjugate:
* $ 
\mathbf{q}^\dagger 
= \begin{bmatrix} q_s \\ -\mathbf{q}_v \end{bmatrix}
$

#### Identity:
  * Unit quaternions: $\|\mathbf{q}\| = 1$
  * $ 
\mathbf{q} \odot \mathbf{q}^\dagger 
= \mathbf{q}^\dagger \odot \mathbf{q}
= \begin{bmatrix} 1 \\ \mathbf{0} \end{bmatrix}
$

#### Rotation action:
* $
Q(\mathbf{q})\mathbf{v} 
= H^T(\mathbf{q} \odot (Hv) \odot \mathbf{q}^\dagger
= H R(\mathbf{q})^T L(\mathbf{q}) H v
$
  * $H = \begin{bmatrix} \mathbf{0} \\ I \end{bmatrix}, \quad Hv = \begin{bmatrix} 0 \\ \mathbf{v} \end{bmatrix}$

#### Rotation correspondence: Passive-Body-to-World (PBTW)
* ${}^W Q^B = Q({}^W q^B)$
* Columns of ${}^W Q^B$ are B-frame axes in W-frame:
  * ${}^W Q^B = \begin{bmatrix} \mathbf{x}^B &\mathbf{y}^B &\mathbf{z}^B \end{bmatrix}$
* Rows of ${}^W Q^B$ are W-frame axes in B-frame:
  * ${}^W Q^B = \begin{bmatrix} (\mathbf{x}^W)^T \\ (\mathbf{y}^W)^T \\ (\mathbf{z}^W)^T \end{bmatrix}$
* Aside: an active (in-frame) rotation uses the same rotation-action definition as PBTW: 
  * $v^W_\text{after} = Q({}^\text{after}\mathbf{q}^\text{before}) v^W_\text{before}$
 
#### Kinematics: Body-frame angular rate
* $ 
{}^W \dot{\mathbf{q}}^B 
= \frac12 {}^W \mathbf{q}^B \odot \begin{bmatrix} 0 \\ \mathbf{\omega}^B \end{bmatrix} 
$

#### Composition:
* $
v^C = {}^C Q^B {}^B Q^A v^A 
= Q({}^C \mathbf{q}^B) Q({}^B \mathbf{q}^A) \mathbf{v}^A 
= Q({}^C \mathbf{q}^B \odot {}^B \mathbf{q}^A) v^A
$

#### Double Cover:
  * $Q(\mathbf{q}) = Q(-\mathbf{q})$
  * Resolution: $\mathbf{q} = \begin{cases} \mathbf{q} & q_s \geq 0 \\ -\mathbf{q} &q_s < 0 \end{cases}$


# Frame Notation
#### Attitude `q_<destination>_<source>`
  * ${}^B \mathbf{q}^A = $ `q_B_A` is the quaternion that represents the PBTW transformation from frame A to frame B.
  * example: ${}^\text{ECI} \mathbf{q}^\text{B} = $ `q_ECI_B` the body-to-inertial state quaternion.

#### Static vectors `u_<expressed in>`
  * $\mathbf{u}^A = $ `u_A` is the vector $\mathbf{u}$ expressed in frame A.
  * example: $\mathbf{\mathrm{magfield}}^B = $ `magfield_B` is the magnetic field vector expressed in the body frame.
 
#### Derivative vectors `v_<of frame>_<with respect to>_<expressed in>`
  * $\mathbf{v}^{C} =$ `v_A_B_C` is the vector $\mathbf{v}$ of frame A, with respect to frame B, expressed in frame C.
  * example: $\mathbf{\omega}^{B} =$ `w_B_ECI_B` the angular velocity of the body frame, with respect to the inertial frame, expressed in the body frame.

# Software Implementations
## Eigen / razzgnc
The following are verified in [`razzgnc/sim/test/test-quaternion-convention.cpp`](https://gitlab.com/albedospace/gnc/razzgnc/razzgnc/-/blob/27df228368dd71f2b1e5e5e7b82eb801b4d677fe/lib/razzgnc/sim/test/test-quaternion-convention.cpp)
#### Scalar/Vector order
  * Scalar first in function interface
    * `auto q = Quaterniond(q_s, q_x, q_y, q_z)`
  * Scalar last in storage
    * `q.data == Vector4d(q_x, q_y, q_z, q_s)`
#### Product: Right-handed/Hamilton
```
Quaterniond i(0, 1, 0, 0);
Quaterniond j(0, 0, 1, 0);
auto ij = (i * j).value();
ij.w() == 0.0; // true
ij.vec() == Vector3d::UnitZ(); // true
```
#### Rotation action: Passive-Body-to-World
  * `v_W = q_W_B * v_B`
  * `v_W = q_W_B.toRotationMatrix() * v_B`
#### Composition: left-to-right
  * `q_C_A = q_C_B * q_B_A`
#### Conjugate
  * `q.inverse()*q == q*q.inverse() == q::Identity()`
#### Normalize
  * `q.normalize()`
  * `p = q.normalized()`

# Literature Comparison
* Matches Sommer.
* Does not match Markley & Crassidis.
* Matches Brian Jackson's Planning with Attitude, though he calls PBTW "Active".
* Matches Sola, though Sola does not distinguish PBTW and "Active".
