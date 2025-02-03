19a20
>         ## いじる必要のないパラメータ
208c209
<     prcp_obs = 5 / (24 * 60 * 60) # kg/m2/s
---
>     prcp_obs = 2 / (24*60*60) # kg/m^2/s
211c212,214
<     pressure=900
---
>     pbtm = 70000
>     ptop = 40000
>     # pbtm-ptop の値しかつかっていないので、cloud_thickness [Pa] とかでもいい
213,216c216,219
<     rho_air = pressure * 100 / (R_d * 255)
<     u_fall = 1.5
<     snow = prcp_obs / (u_fall * rho_air) / (1 - config["resub_factor"]) * 1000 # g/kg
<     delta_snow = rayleigh_step(alpha_final, q_final, snow, delta_final/1000)    
---
>     grav = 9.80665 # m/s^2
>     airmass = (pbtm-ptop)/grav # kg/m^2
> #    u_fall = 1.5 # m/s # もういらないかも
>     dq_snow = prcp_obs/(1-config["resub_factor"])/airmass * 1000 # dq/dt : g/kg/s
218c221,231
<     print(q_final,snow, delta_final, delta_snow*1000)
---
> #    print("init delta = ",delta_final)
>     delta_cloud_vapor = delta_final/1000
>     duration=1
>     delta_snow=0
>     for loop in range(duration*24*60*60):
>         hoge = rayleigh_step(alpha_final, q_final, dq_snow, delta_cloud_vapor)
>         delta_cloud_vapor = ( hoge*(q_final-dq_snow) + delta_final/1000*dq_snow )/q_final
>         delta_snow = delta_snow + ( (hoge + 1) * alpha_final - 1 )/(duration*24*60*60)
> #        print (loop, ( (hoge + 1) * alpha_final - 1 )*1000)
>         
> #     print(q_final,snow, delta_final, delta_snow*1000,R_d)
220,222d232
<     # delta_snow = (rayleigh_results["delta"][-1] / 1000 + 1) * alpha_final - 1 
<     # snow = rayleigh_results["q"][-2] - rayleigh_results["q"][-1] # dummy
< 
228c238
<             snow, 
---
>             dq_snow*duration*24*60*60, # 単位時間あたりではなく、総量を参照すべきと思われる
