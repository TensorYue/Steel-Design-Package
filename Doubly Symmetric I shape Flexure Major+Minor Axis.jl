### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 2d2ea920-b70b-11ee-06fa-8348cf732416
begin
	
	# Material & Total Geometry

	E=29000 # ksi Elastic Modulus
	G=11200 # ksi Shear Modulus
	F_y=50 # ksi Yield Strength
	L_x=18*12 # in
	L_y=18*12 # in
	L_z=18*12 # in
	L_cx=18*12 # in
	L_cy=18*12 # in
	L_cz=18*12 # in
	
	# Cross Section Geometry	
	
	A=72.5 # in^2 Area
	d=36.7 # in Total Depth
	
	t_w=0.8 # in Web Thickness
	
	b_f=16.5 # in Flange Width
	t_f=1.35 # in Flange Thickness

	k=2.3 # in Corner, h=d-2*k for web depth
	h=d-2*k
	
	γ=247/100 # klf Unit Weight over Length

	λ_f= 6.11 # Flange Ratio
	λ_w=40.1 # Web Ratio, can also be determined as (d-2*k)/t_w

	I_x=16700
	I_y=1010

	S_x=913
	S_y=123

	r_x=15.2
	r_y=3.74

	Z_x=1030
	Z_y=190

	r_ts=4.42
	h_0=35.4

	J=34.7
	C_w=316000
end

# ╔═╡ 20961fe2-8ec9-42f7-afb9-87d0504d20d9
# All current ASTM A6 W, S, M, C and MC shapes except W21x48, W14x99, W14x90, W12x65, W10x12, W8x31, W8x10, W6x15, W6x9, W6x8.5, and M4x6 have compact flanges for F_y=50ksi; all current ASTM A6 W, S, M, HP, C, and MC shapes have compact webs at F_y<=70 ksi.

# ╔═╡ 92777bbd-bfc0-44d3-b38f-51678d95f272
47929.91445765297/1.67/12

# ╔═╡ ca495aa0-d5ad-4df3-bc07-e98021ce43e6
# Remember to convert unit (in to ft)

begin

# Major Axis
	
	# Check Compactness
	λ_pf = 0.38*sqrt(E/F_y) # B4.1b
	λ_rf = sqrt(E/F_y) # B4.1b
	λ_pw = 3.76*sqrt(E/F_y) # B4.1b
	λ_rw = 5.7*sqrt(E/F_y) # B4.1b

	c = 1 # F2-8a

	L_b = max(L_x,L_y) # F2 L_b is the length between points that are either braced against lateral displacement of the compression flange or braced against twist of the cross section

	M_p = F_y*Z_x # M_n_Y Yield Resistance F2-1
	C_b = 1 # Lateral Torsional Modification Factor F1. 

	if λ_w <= λ_pw # Compact Web B4.1b

		
		L_p = 1.76*r_y*sqrt(E/F_y) # F2-5
		L_r = 1.95*r_ts*E/0.7/F_y*sqrt(J*c/S_x/h_0 + sqrt((J*c/S_x/h_0)^2+6.76*(0.7*F_y/E)^2)) # F2-6
		
		# Lateral Tortional Buckling
		if  L_b <= L_p
			M_n_LTB = M_p
		elseif L_b <= L_r
			M_n_LTB = C_b*(M_p - (M_p-0.7*F_y*S_x)*(L_b-L_p)/(L_r-L_p)) # F2-2
		else
			F_cr = C_b*pi^2*E/(L_b/r_ts)^2*sqrt(1+0.078*J*c/S_x/h_0*(L_b/r_ts)^2) # F2-4
			M_n_LTB = F_cr*S_x # F2-3
		end

		if λ_f <= λ_pf # Compact Flange B4.1b
			M_n = min(M_p, M_n_LTB) #F2

		# Flange Local Buckling
		elseif λ_f <= λ_rf # Non-Compact Flange
			M_n = min(M_p-(M_p-0.7*F_y*S_x)*(λ_f-λ_pf)/(λ_rf-λ_pf), M_n_LTB, M_p) #F3-1
			
		else # Slender Flange 
			k_c = 4/sqrt(h/t_w)
			M_n =  min(0.9*E*k_c*S_x/λ_f^2, M_n_LTB, M_p) #F3-2
		end
		
	elseif λ_w <= λ_rw # Non-Compact Web B4.1b

		I_yc = 1/12*t_f*b_f^3 # F4.2(c)(2) I_yc is the moment of inertial of the compression flange about the y-axis
		S_xc = (1/12*b_f*t_f^3+b_f*t_f*(0.5*d-0.5*t_f)^2)/(d/2) # F4.1 S_xc is the elastic section modulus referred to compression flange
		h_c = h # F4.2(c)(6) h_c is twicee the distance from the centroid to the following: the inside face of the compression flange less the fillet or corner radius, for rolled shapes; the nearest line of fasteners at the compression flange or the inside faces of the compression flange when welds are used, for built-up sections
		M_yc = F_y*S_xc

		# F4.2(c)(6)
		if I_yc/I_y <= 0.23
			R_pc = 1 # F4-10
		else
			if h_c/t_w <= λ_pw # h_c/t_w
				R_pc = M_p/M_yc # F4-9a
			else
				R_pc = min(M_p/M_yc - (M_p/M_yc-1)*(h_c/t_w-λ_pw)/(λ_rw-λ_pw), M_p/M_yc) # F4-9b
			end
		end

##################
		
		M_n_CFY = R_pc*M_yc # F4-1

		a_w = h_c*t_w/b_f/t_f
		r_t = b_f/sqrt(12*(1+1/6*a_w))

		F_L=0.7*F_y # F4-6a
		
		L_p = 1.1*r_t*sqrt(E/F_y)
		L_r = 1.95*r_t*E/F_L*sqrt(J/S_xc/h_0 + sqrt((J/S_xc/h_0)^2+6.76*(F_L/E)^2))

		if L_b <= L_P
			M_n_LTB = M_p
		elseif L_b <= L_r
			M_n_LTB = C_b*(R_pc*M_yc - (R_pc*M_yc-F_L*S_xc)*(L_b-L_p)/(L_r-L_p))
		else
			F_cr = C_b*pi^2*E/(L_b/r_t)^2*sqrt(1+0.078*J/S_xc/h_0*(L_b/r_t)^2)
			M_n_LTB = F_cr*S_xc
		end

		if λ_f <= λ_pf # Compact Flange B4.1b
			M_n = min(M_n_LTB, M_n_CFY, M_p) # F4

		# Flange Local Buckling
		elseif λ_f <= λ_rf # Non-Compact Flange
			M_n = min(R_pc*M_yc-(R_pc*M_yc-F_L*S_xc)*(λ_f-λ_pf)/(λ_rf-λ_pf), M_n_LTB, M_n_CFY, M_p) # F4-13
			
		else # Slender Flange 
			# F4.3(c)
			if 4/sqrt(h/t_w)<=0.35
				k_c = 0.35
			else
				k_c = min(0.76, 4/sqrt(h/t_w))
			end
			M_n =  min(0.9*E*k_c*S_x/(b_f/2/t_f)^2, M_n_LTB, M_n_CFY, M_p) # F4-14
		end

		# Doubly Symmetric Shape, no tension flange yielding

####################
	
	else # Slender Web B4.1b

		I_yc = 1/12*t_f*b_f^3 # F4.2(c)(2) I_yc is the moment of inertial of the compression flange about the y-axis
		S_xc = (1/12*b_f*t_f^3+b_f*t_f*(0.5*d-0.5*t_f)^2)/(d/2) # F4.1 S_xc is the elastic section modulus referred to compression flange
		h_c = h # F4.2(c)(6) h_c is twicee the distance from the centroid to the following: the inside face of the compression flange less the fillet or corner radius, for rolled shapes; the nearest line of fasteners at the compression flange or the inside faces of the compression flange when welds are used, for built-up sections

		a_w = min(h_c*t_w/b_f/t_f, 10) # F5.2(c) a_w is defined by Eq F4-12, but shall not exceed 10
		r_t = b_f/sqrt(12*(1+1/6*a_w))

		L_p = 1.1*r_t*sqrt(E/F_y)
		L_r = pi*r_t*sqrt(E/0.7/F_y)

		R_pg = min(1-a_w/(1200+300*a_w)*(h_c/t_w-5.7*sqrt(E/F_y)), 1)

		M_n_CFY = R_pg*F_y*S_xc # F5-1

		if L_b <= L_p
			F_cr = F_y
		elseif L_b <= L_r
			F_cr = min(C_b*(F_y - (0.3*F_y)*(L_b-L_p)/(L_r-L_p)), F_y)
		else
			F_cr = min(C_b*pi^2*E/(L_b/r_t)^2, F_y)
		end
		
		M_n_LTB = R_pg*F_cr*S_xc # F5-2

		if λ_f <= λ_pf # Compact Flange B4.1b
			M_n = min(M_n_LTB, M_n_CFY, M_p) # F5

		# Flange Local Buckling
		elseif λ_f <= λ_rf # Non-Compact Flange
			M_n = min(R_pg*(F_y-0.3*F_y*(λ_f-λ_pf)/(λ_rf-λ_pf))*S_xc, M_n_LTB, M_n_CFY, M_p) # F5-8
			
		else # Slender Flange 
			if 4/sqrt(h/t_w)<=0.35
				k_c = 0.35
			else
				k_c = min(0.76, 4/sqrt(h/t_w))
			end
			M_n =  min(R_pg*0.9*E*k_c*S_x/(b_f/2/t_f)^2, M_n_LTB, M_n_CFY, M_p) # F5-9
		end

		# Doubly Symmetric Shape, no tension flange yielding
	end

	
	print("remember to convert unit")
	M_n # Unit in
	M_n/12 # Unit ft

# Minor Axis

	m_p = min(F_y*Z_y, 1.6*F_y*S_y)

	if λ_f <= λ_pf # Compact Flange B4.1b
		m_n = m_p #F2

	# Flange Local Buckling
	elseif λ_f <= λ_rf # Non-Compact Flange
		m_n = min(m_p-(m_p-0.7*F_y*S_y)*(λ_f-λ_pf)/(λ_rf-λ_pf), m_p) #F3-1
			
	else # Slender Flange 
		m_n =  min(0.69*E*S_y/λ_f^2, m_p) #F3-2
	end

	m_n # Unit in
	m_n/12 # Unit ft
	
end
		

# ╔═╡ 15205757-4c57-4dc8-90dd-58ead78f5a31
# Minor Axis!!!

begin
	

# ╔═╡ Cell order:
# ╠═2d2ea920-b70b-11ee-06fa-8348cf732416
# ╠═20961fe2-8ec9-42f7-afb9-87d0504d20d9
# ╠═92777bbd-bfc0-44d3-b38f-51678d95f272
# ╠═ca495aa0-d5ad-4df3-bc07-e98021ce43e6
# ╠═15205757-4c57-4dc8-90dd-58ead78f5a31
