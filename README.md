# Monte Carlo Path Tracing 

## How it's Made:
This project uses the Julia language and a couple Julia libraries. 

## Build and Run
1. ```cd``` into the ```RayTracer``` directory
2. ```julia --project``` to open the Julia CLI
3. For development: ```using revise``` (optional)
4. ```using RayTracer```
5. ```RayTracer.main(scene, camera, width, sample_couunt, outfile)```

## Rendering Algorithms
This began using the simple Whitted style algorithm where rays are traced from the viewpoint into the scene and rays are recursively traced when a reflective or refractive surface is interacted with. This is a quick way to get something up and running to produce decent looking images, but it lacks the fidelity to create certain lighting effects such as glossy reflection, soft shadows, and translucency. To address these shortcomings, I implemented distributed ray tracing methods following Cook et. al., 1984. In a nutshell, multiple rays are traced per pixel, and rather than tracing the ideal reflection or refaction direction upon those types of interactions, a probability distribution is sampled from to randomly generate reflection and refraction rays around the ideal directions. The only thing missing from distributed ray tracing is indirect lighting interactions where surfaces are illuminated by light that reflects off one or more surfaces before hitting the surface of interest. To implement indirect lighting, an algorithm known as Monte Carlo Path Tracing is used. Multiple rays are still sent from each pixel and rather than stopping when a diffuse surface is reached, the hemisphere over the intersected point is sampled according to some probability distribution, a new ray is generated using this sampling and then recursively traced. This process repeats until a ray intersects a light source, at which point the final lighting value is known for that path sample. 

## Optimizations:
Due to the randomized nature of Monte Carlo Path tracing, images are often noisy when the sample count is low. The easy solution is to generate more rays per pixel, which reduces noise at the cost of making renders take significantly longer. There are multiple approaches to noise reduction. The first being importance sampling: generate samples in the hemisphere where a larger lighting contribution will be made. To this end I implemented the standard cosine-weighted hemisphere sampling which generates proportionally fewer rays at the edges of the hemisphere. compared with unifomly sampling the hemisphere, the expensive samples are being used more efficently. The second appraoch to noise reduction that I implemented is known as next event estimation. This is a method to compute direct lighting values for surface interactions by combining samples from multiple distributions. Direct lighting can be found in two ways: a ray generated randomly in the hemisphere travels directly to a light source on the next bounce, or a point on the light source is sampled directly and a shadow ray is traced from the intersection point to the sampled light source point to determine if the surface point is illuminated by the sampled light source. Combining next event estimation with importance sampling produces a much clearer image using less samples than if these techniques were not used. 

## Key Takeaways:
Designing a robust API for materials, lighting interactions, and the many PDFs that need to be sampled was quite the challenge. At first all lighting was done in a single ```determine_color``` method which was not at all extensible to distributed ray tracing or path tracing. Julia has a neat feature called multiple dispatch where instead of overloading a method with different arguments, the method is overloaded with different object types. Using this feature I experimented with a set of ```determine_color``` methods, each of which dispatched on a different material type. To transition to path tracing, I created a set of ```sample, eval, and pdf``` methods that each dispatched on material type so no matter the type of material being intersected, all that needed to be done is generate the next bounce direction using ```sample``` and compute the current surface's "lighting contribution" using ```eval and pdf``` methods. 

## Future Work
* Renders seem to be slightly too dark, fix this bug in lighting calculations.
* Extend the materials API to handle dielectric materials (again)
* Port this thing into C++ to render things faster (see my Light-Transport repo)
