
# The Great Finite Volume Library

The FVLib code is a library of advanced computational algorithms and numerical methods to solve partial differential equations within the **finite volume philosophy**. The project aims to deliver **high-accurate**, **high-performance**, and **high-efficient** simulations of a wide range of physics and mechanics problems in relevant aerospace, manufacturing, environmental, and biomedical applications.

The FVLib code is the result of years of dedication and passion for applied mathematics and scientific computing, and the ambition of **pushing the limits of numerical simulation even further** in serving the scientific knowledge and technological innovation.

---

## Modern object-oriented Fortran (2003/2008 standards)

The FVLib code is programmed in modern Fortran (2003/2008 standards) with an object-oriented paradigm for better **code reuse, maintenance, and readability**. Its architecture is organised in three levels:

- **Core level**: including linear algebra algorithms, input/output interfaces, mesh and field handlers, sparse matrix structures, etc.

- **Applications level**: including specific model solvers, specific boundary and interface conditions, and pre-processing and post-processing tools.

- **Cases level**: including geometry files, mesh generation scripts, models and schemes parameters files, and running and post-processing scripts.

<br>

<div align="center">
  <table>
    <tr>
      <td align="center"><img src="images/applications_level.png" width="100%"></td>
      <td align="center"><img src="images/cases_level.png" width="100%"></td>
    </tr>
    <tr>
      <td align="center">Snapshot of applications level.</td>
      <td align="center">Snapshot of cases setup.</td>
    </tr>
  </table>
</div>

---

## Highly accurate schemes in space and time

The discretisation methods implemented in the FVLib code are highly accurate in space and time, effectively achieving up to the eighth-order of convergence. Comprehensive benchmarking proves that high-order accurate schemes benefit from a **better trade-off between accuracy and efficiency** than the counterpart lower-order accurate ones. This property can be exploited in different ways:

- **Improved accuracy**: for the same discrete geometrical representation level (number of degrees of freedom), high-order accurate schemes provide **significantly more accurate solutions** than those obtained with the traditional first- and second-order accurate schemes.

- **Enhanced performance**: for the same approximate solution accuracy level, high-order accurate schemes provide **significantly more efficient computations** (execution time) than those of the traditional first- and second-order accurate schemes.

- **Resource-use efficiency**: for the same approximate solution accuracy level, high-order accurate schemes consume **significantly less resources** (power and memory) than those of the traditional first- and second-order accurate schemes.

<br>

<div align="center">
  <table>
    <tr>
      <td align="center"><img src="images/error_vs_dof.png" width="70%"></td>
      <td align="center"><img src="images/error_vs_time.png" width="70%"></td>
    </tr>
    <tr>
      <td align="center">Error <em>versus</em> number of unknowns.</td>
      <td align="center">Error <em>versus</em> execution time.</td>
    </tr>
  </table>
</div>

---

## Unstructured meshes for complex geometries

Complex geometries arise in many real-world problems of physics and engineering applications, for which domain fitted unstructured meshes are still the preferred approach for its **flexibility and robustness**, especially in 3D.

- **General element shapes**: the high-order accurate discretisation methods implemented in the FVLib code can handle **2D and 3D unstructured meshes with general element shapes** for the most demanding problems in intricate geometries.

- **Piecewise linear meshes**: the high-order accurate discretisation methods implemented in the FVLib code preserve the optimal **high-orders of convergence with the standard linear piecewise elements**, on arbitrary curved geometries, overcoming the cumbersomeness of generating and dealing with curved meshes as the traditional approaches.

<br>

<div align="center">
  <table>
    <tr>
      <td align="center"><img src="images/curved_domain.png" width="70%"></td>
      <td align="center"><img src="images/unstructured_mesh.png" width="70%"></td>
    </tr>
    <tr>
      <td align="center">Intricate curved geometry.</td>
      <td align="center">Unstructured piecewise linear mesh.</td>
    </tr>
  </table>
</div>

- **Gmsh compatible interface**: the FVLib code provides a compatibility interface with [Gmsh](https://gmsh.info/), an open source **3D mesh generator with built-in pre- and post-processing facilities**, designed to provide a fast, light and user-friendly meshing tool with parametric input and flexible visualisation capabilities.

<br>

<div align="center">
  <table>
    <tr>
      <td align="center"><img src="images/gmsh1.png" width="100%"></td>
      <td align="center"><img src="images/gmsh2.png" width="100%"></td>
    </tr>
    <tr>
      <td align="center">Snapshot of Gmsh with geometry.</td>
      <td align="center">Snapshot of Gmsh with mesh.</td>
    </tr>
  </table>
</div>

<br>

<div align="center">
  <table>
    <tr>
      <td align="center"><img src="images/gmsh3.png" width="100%"></td>
      <td align="center"><img src="images/gmsh4.png" width="100%"></td>
    </tr>
    <tr>
      <td align="center">Snapshot of Gmsh with solution.</td>
      <td align="center">Snapshot of Gmsh with post-processing.</td>
    </tr>
  </table>
</div>

---

## Parallel computing for HPC environments

The FVLib code multi-processing capabilities are implemented based on the **shared-memory and message-passing parallel programming models** through the OpenMP and OpenMPI application programming interfaces. It allows the development and deployment of portable and scalable large-scale parallel applications that take advantage of modern HPC systems.

---

## Contributing

The FVLib code is not currently an open-source project. However, anyone willing to contribute to the project and/or making use of its capabilities on a collaboration basis is welcome. If you are interested in using it or need further details, you can always [contact me](mailto:rcosta@dep.uminho.pt).
