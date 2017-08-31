readcol, 'kins_meanspread.txt', mid_a, mid_r, top_a, top_r, bottom_a, bottom_r, midtop_a, midtop_r, mibottom_a, midbottom_r, time, $
       f='D,D,D,D,D,D,D,D,D,D,A'

       utbasedata = anytim(time[0])
       t = anytim(time) - utbasedata
       rsun = 6.95508e8

       mid_km = mid_r * rsun / 1000.0d
       top_km = top_r * rsun / 1000.0d
       bottom_km = bottom_r * rsun / 1000.0d
       midtop_km = midtop_r * rsun / 1000.0d
       midbottom_km = midbottom_r * rsun / 1000.0d

       mid_vel = deriv(t, mid_km)
       top_vel = deriv(t, top_km)
       bottom_vel = deriv(t, bottom_km)
       midtop_vel = deriv(t, midtop_km)
       midbottom_vel = deriv(t, midbottom_km)

       mid_accel = deriv(t, mid_vel) * 1000.0d ; to be in Metres
       top_accel = deriv(t, top_vel) * 1000.0d
       bottom_accel = deriv(t, bottom_vel) * 1000.0d
       midtop_accel = deriv(t, midtop_vel) * 1000.0d
       midbottom_accel = deriv(t, midbottom_vel) * 1000.0d
     
end       