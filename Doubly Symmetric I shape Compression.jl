### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 68ae7da0-b6f6-11ee-0741-95699353f894
# Doubly Symmetric I Shape Wale Design

# ╔═╡ af5fd033-07ef-422f-9187-4369561d6bb3
# Input Geometry
begin
	
	# Material & Total Geometry

	E=29000 # ksi Elastic Modulus
	G=11200 # ksi Shear Modulus
	F_y=50 # ksi Yield Strength
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

# ╔═╡ 0d744723-d947-49d0-8eab-7da17d00d163
# Compact Check
begin

	λ_c_rf = 0.56*sqrt(E/F_y) # B4.1(a)
	
	λ_c_rw = 1.49*sqrt(E/F_y) # B4.1(a)

	λ_lr = max(L_cx/r_x, L_cy/r_y) # E3-2, E3-3
	
	# Determine F_e (Flexural Buckling / Tortional Buckling)

	# C_w = I_y*h_0^2/4 # E4(c)
	F_e = min(pi^2*E/(λ_lr)^2, (pi^2*E*C_w/L_cz^2+G*J)/(I_x+I_y)) #  E3-4, E4-2
	
	# Determine F_cr (From F_e)
	
	if λ_lr <= 4.71*sqrt(E/F_y) # E3(a)
			F_cr = 0.658^(F_y/F_e)*F_y # E3-2
	else # E3(b)
			F_cr = 0.877*F_e # E3-3
	end
	
	# Determine A_ge (Non-Slender / Slender)

	if λ_f <= λ_c_rf && λ_w <= λ_c_rw # E3, E4, E7
		A_ge = A
	else

		h = λ_w*t_w
		
		A_flange = 2*b_f*t_f
		A_web = h*t_w
		A_res = A - A_flange - A_web # Area except Web and Flange
		
		# Flange

		if λ_f <= λ_c_rf*sqrt(F_y/F_cr) # E7(a)
			A_flange = A_flange # E7-2
		else # E7(b)
			c_1_f = 0.22 # E7.1
			c_2_f = (1-sqrt(1-4*c_1_f))/2/c_1_f # E7-4
			F_el_f = (c_2_f*λ_c_rf/λ_f)^2*F_y # E7-5

			A_flange = A_flange * (1-c_1_f*sqrt(F_el_f/F_cr))*sqrt(F_el_f/F_cr) # E7-3
		end

		# Web

		if λ_w <= λ_c_rw*sqrt(F_y/F_cr) # E7(a)
			A_web = A_web # E7-2
		else
			c_1_w = 0.18 # E7.1
			c_2_w = (1-sqrt(1-4*c_1_w))/2/c_1_w # E7-4
			F_el_w = (c_2_w*λ_c_rw/λ_w)^2*F_y # E7-5

			A_web = A_web * (1-c_1_w*sqrt(F_el_w/F_cr))*sqrt(F_el_w/F_cr) # E7-3
		end

		A_ge = A_res + A_flange + A_web
		
	end

	P_n = F_cr*A_ge # E3-1, E4-1, E7-1
	print(A_ge)
end

# ╔═╡ f4c16b6c-d8c0-41e7-818b-097e77c01be5
begin
	a=1
	b=2
	if a<=1 && b<=1
		print(1)
	else
		print(0)
	end
end

# ╔═╡ 3741e2ae-89af-44b0-8aef-cdc1a08bb973
begin
	c=2
	c=c+1
end

# ╔═╡ Cell order:
# ╠═68ae7da0-b6f6-11ee-0741-95699353f894
# ╠═af5fd033-07ef-422f-9187-4369561d6bb3
# ╠═0d744723-d947-49d0-8eab-7da17d00d163
# ╠═f4c16b6c-d8c0-41e7-818b-097e77c01be5
# ╠═3741e2ae-89af-44b0-8aef-cdc1a08bb973
