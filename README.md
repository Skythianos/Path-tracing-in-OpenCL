# Path-tracing-in-OpenCL
It is significantly easier and more cost effective to render photorealistic images of 3D models created by 3D modeler/designer programs (architectural, automotive etc.), than building them up, just to see how they would look in real life.
But, physically based rendering is not a trivial task, and it is pretty computing-intensive. Algorithms like this are called global illumination algorithms. In this project, I will present the most basic algorithm among them, called path tracing, in OpenCL.
To get a physically correct image, it raises some questions. How to model diffuse
surfaces? How do they reflect the rays? How does an ideal mirror look? How refractive
materials can be modeled? What are the Fresnel equations?
The task is very computing-intensive, so it would be good to speed up the
calculations. I use two things in this project: brain and brute force. On the one hand, with
the bounding boxes of the objects, we can determine very efficiently (sometimes with
only one intersection test), if the ray will intersect an object with even millions of
triangles. On the other hand, GPUs nowadays outperforms CPUs in massively parallel
problems, so it is not a surprise, that I will use the GPU to do the calculations.
In this project, I make a suggestion to decompose the problem to independent
threads, I show the speed up by using the GPU instead of the CPU to render the images,
I show a working solution that how can the CPU and the GPU work cooperatively in
OpenCL, and I will answer the questions above as well.
