Path Tracer
======================

University of Pennsylvania, CIS 561: Advanced Computer Graphics, Course project
------------


Features :
------------

- Normal Map

- Diffuse are light, point light, spot light, infinite light(Environment Map)

- Different materials (Metal, ChromMirror, glass, plastic, etc.)

- Acceleration : BVH / Kd-Tree

- Lens-based Camera

- Implicit Surface Rendering

- Photon Mapper (may take super longer time, but render merges quite well!)

- Volume Rendering



Rendered Images
------------

Volume Rendering :

Original Photon Mapper render

![](./renders/hw11_100_400SamplesPM.png)

thin & thick white fog volume render 

![](./renders/hw11_rendered_images7_160000.png)  ![](./renders/hw11_rendered_images7_360000.png)  


final volume render we get

![](./renders/hw11_100_400SamplesPM_thin_fog.jpg)  ![](./renders/hw11_100_400SamplesPM_heavy_fog.jpg)  





Photon Mapper : Rainbow Box

![](./renders/hw10_61_400SamplePM.png)




Implicit Surface Rendering using Ray-marching

![](./renders/hw9_implicit_surface_62_test8.png)

![](./renders/hw9_implicit_surface_66.png)

![](./renders/hw9_implicit_surface_68_6.png)



Lens-based Camera with Lens Radius = 1.0 , 5.0, 10.0

![](./renders/hw9_thin_len_60_lensR_1.0_f_30.png)  

![](./renders/hw9_thin_len_60_lensR_5.0_f_30.png)

![](./renders/hw9_thin_len_60_lensR_10.0_f_30.png)



Complex model rendered within just few minutes (general <10min depends on your computer)with BVH/Kd-Tree and more the 1 hour at least without it.

900 ssp MIS Glass dragon (around 3min rendering time) / 900 ssp MIS Happy Buddha (around 4min rendering time)

![](./renders/hw8_accel_47_900SampleMIS_8recur.png)  ![](./renders/hw8_accel_48_900SampleMIS_8recur.png) 







victory pose (Mirror Material)

![](./renders/53_900MIS_5recur.png) 



Glass cup and ball (Environment Map)

![](./renders/hw10_69_900SampleMIS.png) 


AirJordan Jumpman (Glass Material + Environment Map)

![](./renders/63_4.png)

(This one actually a little bit light / uv problem, but I think it's a cool mistake)

![](./renders/63.png)



Mountains (Texture + Normal Map + Environment Map)

![](./renders/66.png)



Wahoo series

![](./renders/55_2.png)  ![](./renders/55.png)


Wahoo's War(several different color point lights)

![](./renders/hw10_66_900SampleMIS.png)


Evil Wahoo..(OK, it's just a light problem, too dark)

![](./renders/54_uvProblem.png)



Normal Map Comparision (left back wall with normal map, right without)

![](./renders/57_100MISWithNormalMap.png)  ![](./renders/57_100MIS.png)



TwoLights Cornell Box

![](./renders/52_5Recursion.png)



Veach scene with Naive and direct lighting integrator

![](./renders/5_900Naive_Veach.png)   ![](./renders/5_900Direct_Veach.png) 


Glass Ball scene

![](./renders/64_eta_1.5.png.png)




