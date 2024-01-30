### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 17f8d8de-6fd2-423b-a128-6d31bc3d7d8c
function Round_HSS_Design(D,t,L,F_y,P,M)
	D_i=D-2*t
	A=pi/4*(D^2-D_i^2)
	S=pi/32*(D^4-D_i^4)/D
	Z=1/6*(D^3-D_i^3)
	I=pi/64*(D^4-D_i^4)
	r=sqrt(I/A)
	
	E=29000
	λ=D/t
	
	# Compression
	
	if λ<0.11*E/F_y
		A_ge=A
	else
		A_ge=(0.038*E/F_y/λ+2/3)*A
	end

	L_c=L*12
	F_e=pi^2*E/(L_c/r)^2

	if L_c/r<=4.71*sqrt(E/F_y)
		F_cr=0.658^(F_y/F_e)*F_y
	else
		F_cr=0.877*F_e
	end
	return([A_ge])
	P_n=F_cr*A_ge

	# Flexure
	if λ<0.07*E/F_y
		M_n=F_y*Z/12
	elseif λ<0.31*E/F_y
		M_n=min(F_y*Z,(0.021*E/λ+F_y)*S)/12
	else
		M_n=min(F_y*Z,0.33*E/λ*S)/12
	end

	P_c=P_n/1.67
	M_c=M_n/1.67
	

	# print(P_c)
	# print(" ")
	# print(M_c)

	if P/P_c>=0.2
	#	return P/P_c+8/9*M/M_c
	else
	#	return P/P_c/2+M/M_c
	end
end

# ╔═╡ bdabe2a1-e6d4-4928-95c4-3f40e1833674
Round_HSS_Design(36,0.4,89.6,45,763.2,306.56)

# ╔═╡ 52e2eef1-14a4-40c1-9fff-b9f3a73d8ae0
begin
	Store=[]
	for i=1:10
		append!(Store, Round_HSS_Design(36,i*0.1,89.6,45,763.2,306.56))
	end
end

# ╔═╡ 9743c0ee-c184-42d3-9849-b0e7f7bd9315
Store

# ╔═╡ Cell order:
# ╠═bdabe2a1-e6d4-4928-95c4-3f40e1833674
# ╠═52e2eef1-14a4-40c1-9fff-b9f3a73d8ae0
# ╠═9743c0ee-c184-42d3-9849-b0e7f7bd9315
# ╠═17f8d8de-6fd2-423b-a128-6d31bc3d7d8c
