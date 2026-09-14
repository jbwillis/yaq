import marimo

__generated_with = "0.24.1"
app = marimo.App()


@app.cell
def _():
    import marimo as mo

    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # YAQ
    Quaternions are a powerful tool for encoding 3D rotations. But, in their definition and use there are many assumptions that lead to competing conventions. This writeup does not enumerate all of the conventions, but it does attempt to be explicit about the specific assumptions made in establishing a convention.

    ## Hyper-Complex Numbers
    There are many ways to introduce quaternions, but I think the best could be the way they were conceived: as hyper-complex numbers.
    A quaternion $\mathbf{q}$ is a hyper-complex number
    $$
    \mathbf{q} = q_s + q_x\mathbf{i} + q_y\mathbf{j} + q_z \mathbf{k}
    $$
    where the following identity governs multiplication
    $$
    \mathbf{i}^2 = \mathbf{j}^2 = \mathbf{k}^2 = \mathbf{i}\mathbf{j}\mathbf{k} = -1.
    $$

    /// note | The order of multiplication matters.
    We chose $\mathbf{i}\mathbf{j}\mathbf{k} = -1$ as opposed to another ordering. This forms a right handed coordinate system, as will be seen shortly.
    ///

    This definition also leads to the following identities:
    $$
    \begin{aligned}
    \mathbf{i}\mathbf{j} &= -\mathbf{j}\mathbf{i} = \mathbf{k}\\
    \mathbf{j}\mathbf{k} &= -\mathbf{k}\mathbf{j} = \mathbf{i}\\
    \mathbf{k}\mathbf{i} &= -\mathbf{i}\mathbf{k} = \mathbf{j}.
    \end{aligned}
    $$

    ## Vector-form
    Each term in the hyper-complex number can be mapped to an element in a 4-vector:

    $$
    \mathbf{q} = \begin{bmatrix} q_s \\ q_x \\ q_y \\ q_z \end{bmatrix} = \begin{bmatrix} q_s \\ \mathbf{q}_v \end{bmatrix}.
    $$
    where we refer to the complex portion of a quaternion as its _vector part_:
    $$
    \mathbf{q}_v = \begin{bmatrix}  q_x \\ q_y \\ q_z \end{bmatrix}
    $$
    and the real member as its _scalar part_: $q_s$.
    A quaternion with $q_s = 0$ is referred to as a _pure quaternion_:
    $$
    \mathbf{q} = \begin{bmatrix} 0 \\ \mathbf{q}_v \end{bmatrix}.
    $$

    /// note | Scalar part first.
    We chose to place $q_s$ at the _top_ of $\mathbf{q}$ to match the common real-then-imaginary ordering of complex numbers.
    ///

    ## The Quaternion Product
    The first thing we will do with the hyper-complex number definition is investigate its behavior when multiplying the vector-part of quaternions. Let
    $$
    \mathbf{p}_v = p_x\mathbf{i} + p_y\mathbf{j} + p_z \mathbf{k} \\
    \mathbf{q}_v = q_x\mathbf{i} + q_y\mathbf{j} + q_z \mathbf{k}.
    $$
    Multiplying following the hyper-complex number rules above,
    $$
    \begin{aligned}
    \mathbf{p}_v * q_v =
    &  p_x\mathbf{i}(q_x\mathbf{i} + q_y\mathbf{j} + q_z \mathbf{k}) \\
    &+  p_y\mathbf{j}(q_x\mathbf{i} + q_y\mathbf{j} + q_z \mathbf{k}) \\
    &+  p_z\mathbf{k}(q_x\mathbf{i} + q_y\mathbf{j} + q_z \mathbf{k}) \\
    =& p_xq_x\mathbf{i}\mathbf{i} + p_xq_y\mathbf{i}\mathbf{j} + p_xq_z\mathbf{i}\mathbf{k} \\
    &+ p_yq_x\mathbf{j}\mathbf{i} + p_yq_y\mathbf{j}\mathbf{j} + p_yq_z\mathbf{j}\mathbf{k} \\
    &+ p_zq_x\mathbf{k}\mathbf{i} + p_zq_y\mathbf{k}\mathbf{j} + p_zq_z\mathbf{k}\mathbf{k} \\
    =& -p_xq_x + p_xq_y\mathbf{k} - p_xq_z\mathbf{j} \\
    &  -p_yq_x\mathbf{k} - p_yq_y + p_yq_z\mathbf{i} \\
    &  +p_zq_x\mathbf{j} - p_zq_y\mathbf{i} - p_zq_z \\
    =& -p_xq_x - p_yq_y - p_zq_z \\
    & -p_zq_y\mathbf{i} +p_yq_z\mathbf{i} \\
    & +p_zq_x\mathbf{j} -p_xq_z\mathbf{j} \\
    & -p_yq_x\mathbf{k} +p_xq_y\mathbf{k}. \\
    \end{aligned}
    $$
    Now, we switch into a vector form to give some structure to the result:
    $$
    \begin{aligned}
    =&
    \begin{bmatrix}
    -p_xq_x - p_yq_y - p_zq_z \\
    -p_zq_y +p_yq_z \\
    +p_zq_x -p_xq_z \\
    -p_yq_x +p_xq_y
    \end{bmatrix}
    =&
    \begin{bmatrix}
    -p_x &-p_y &-p_z \\
    0 &-p_z &p_y \\
    p_z & 0 & -p_x \\
    -p_y &p_x & 0
    \end{bmatrix}
    \begin{bmatrix}
    q_x \\ q_y \\ q_z
    \end{bmatrix}
    =&
    \begin{bmatrix}
    -\mathbf{p}_v^T \\
    [\mathbf{p}_v]_\times
    \end{bmatrix}
    \mathbf{q}_v.
    \end{aligned}
    $$
    The scalar part of the result is the negative dot product of the two vectors, and the vector part of the result is the cross product of the two vectors!
    Interestingly, the quaternion product predates the terms "vector" and "cross-product."

    With a bit of notational abuse, we'll use the above result to define the full quaternion product as a matrix multiply:
    $$
    \begin{aligned}
    \mathbf{p} \odot \mathbf{q}
    &= (p_s + \mathbf{p}_v) * (q_s + \mathbf{q}_v) \\
    &= p_sq_s + p_s \mathbf{q}_v + q_s \mathbf{p}_v + \mathbf{p}_v * \mathbf{q}_v \\
    &=
    \begin{bmatrix}
    p_s & \mathbf{0}^T \\
    \mathbf{p}_v^T & p_s I
    \end{bmatrix}
    \begin{bmatrix} q_s \\ \mathbf{q}_v \end{bmatrix}
    +
    \begin{bmatrix}
    0 & -\mathbf{p}_v^T \\
    \mathbf{0} & [\mathbf{p}_v]_\times
    \end{bmatrix}
    \begin{bmatrix} q_s \\ \mathbf{q}_v \end{bmatrix}\\
    &=
    \begin{bmatrix}
    p_s & -\mathbf{p}_v^T \\
    \mathbf{p}_v^T & p_s I +[\mathbf{p}_v]_\times
    \end{bmatrix}
    \mathbf{q}
    \\
    &\triangleq L(\mathbf{p})\mathbf{q}.
    \end{aligned}
    $$
    """)
    return


@app.cell
def _(mo):
    # Render an image from a URL
    mo.image(
        src="https://commons.wikimedia.org/wiki/File:William_Rowan_Hamilton_Plaque_-_geograph.org.uk_-_347941.jpg",
        # alt="The plaque commemorating Hamilton's breakthrough",
        width=300,
        rounded=True,
        caption="The plaque commemorating Hamilton's breakthrough",
    )
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
