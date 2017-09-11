import numpy as np

cmaps = dict()

# the information is got from Slicer3D presets.xml
ct_aaa_controls = np.array([-3024,143.556,166.222,214.389,419.736 ,3071])
ct_aaa_colors=np.array([[0,0,0,0],
                       [0.615686,0.356863,0.184314,0],
                       [0.882353,0.603922,0.290196,.0], # 0.686275
                       [1,1,1,.0], #0.696078
                       [1,0.937033,0.954531,0.833333],
                       [0.827451,0.658824,1,0.803922]])
ct_aaa_gradient_opacities = np.array([4,0,1,255,1])

cmaps['ct_aaa'] = {'controls':ct_aaa_controls,'colors':ct_aaa_colors,'gradient_opacity':ct_aaa_gradient_opacities}