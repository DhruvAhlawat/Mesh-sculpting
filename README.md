# Mesh Scuptor software
--- 
can import any file in .OBJ format, and then allows to perform extrude/smoothen/subdivision/revolution operations. 
More examples of meshes [here](a2_report.pdf)

## Creative Example

![image](https://github.com/user-attachments/assets/d99ea62a-9199-4301-b594-d7c39beba785)
Image of a render of a stylized spinning top mesh, created by extruding and revolving faces of a sphere multiple times. 

![image](https://github.com/user-attachments/assets/2040299c-48e0-4077-8400-edbb53a0d5f5)
Top view of same.


# Running

Make sure that [glm](https://github.com/g-truc/glm) and [SDL2](https://www.libsdl.org/) are installed. Ideally, these should be installed by your package manager rather than manually (at least, if you are on Linux or Mac).

Then compile the code using the standard CMake procedure:

- The first time, run `cmake -B build` from the project root to create a `build/` directory and initialize a build system there.
- Then, every time you want to compile the code, run `cmake --build build` (again from the project root). Then the example programs will be created under `build/`.
