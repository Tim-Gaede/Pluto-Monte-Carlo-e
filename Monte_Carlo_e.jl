### A Pluto.jl notebook ###
# v0.11.14

using Markdown
using InteractiveUtils

# ╔═╡ ec55d1a4-f2a4-11ea-21ba-d5d9e9fab3e5
begin
	
using Pkg
using Formatting
#using OffsetArrays

end

# ╔═╡ 30ef477e-f33a-11ea-257f-fdb9ba8bf110
e = exp(1)

# ╔═╡ fa15ffc8-f2a2-11ea-39cd-27af9b2559d0
function MCe_targets(rows::Int, cols::Int) 
	
    hits_by_row_col = zeros(Int8, rows, cols) 
   
    for dart in (1 : rows * cols) 		
		row = rand(1:rows)
		col = rand(1:cols)
        hits_by_row_col[row, col] += 1  		
    end
    
	
    hits_by_row_col 
end

# ╔═╡ 171ee65e-f33c-11ea-012e-3b4e1ba23fc9
hits_by_rc = MCe_targets(10, 10)

# ╔═╡ 97cf78d8-f33a-11ea-17fa-6157a0fad0b1
function hit_cnts_from_targets(hits_by_row_col)
	
	rows = size(hits_by_row_col, 1)
    cols = size(hits_by_row_col, 2)
     
	
	hits_max = Int8(0)
	for row in (1:rows)
		for col in (1:cols)
			if hits_by_row_col[row, col] > hits_max
				hits_max = hits_by_row_col[row, col]
			end
		end
    end
    
    
		 
	#result = OffsetArray(zeros(Int64, hits_max+1), 0:hits_max)
	# Create array of zeros with indicies from 0 to hits_max
	result = zeros(Int64, 0:hits_max) 
 	
	for row in (1:rows)
		for col in (1:cols)
			result[hits_by_row_col[row, col]] += 1
		end
	end
	
	
	result
end

# ╔═╡ 3eaca8a0-f33c-11ea-14f7-258438f5b563
hit_cnts = hit_cnts_from_targets(hits_by_rc)

# ╔═╡ 750ae004-f339-11ea-2fdc-3386c26e35a0
function e_ests_from_hit_cnts(a)
 
	ests = []
	
    ∑ = sum(a)
	
    for i in (0 : length(a)-1)
         push!(ests, ∑ /(a[i] * factorial(i)))
    end
     

    ests
end 

# ╔═╡ a7ea64ae-f33e-11ea-00ba-d7e2dd0003ce
e_ests_from_hit_cnts(hit_cnts)

# ╔═╡ b38f332c-f328-11ea-3ab7-c5ed3b530bbb
# Add more to the array of number of "tets"
# A "tet" is a set (single, pair, trio, quartet, quintet, hextet, heptet, etc.) 
# of random variables (0.0 to 1.0) whose sum exceeded 1.0



function MCe_∑rands!(numTets_by_index, numTets_additional) 
    
    for t = 1 : numTets_additional
        # Sum random variables (0 to 1) until a sum of > 1 is reached
        # Count the number of random variables
        ∑ = 0.0
        cnt = 0
        while ∑ ≤ 1.0
            ∑ += rand()
            cnt += 1
        end

       
        if cnt <= length(numTets_by_index)
            numTets_by_index[cnt] += 1
        else
            # Extend the array
            num_Intermediates = cnt - length(numTets_by_index) - 1
            for newIter in (1 : num_Intermediates)
                push!(numTets_by_index, 0)
            end
            push!(numTets_by_index, 1)
        end        
    end


    numTets_by_index
end

# ╔═╡ d5a282a6-f333-11ea-209f-2d6d48671f3e
function e_est_from_tets(numTets_by_index)
    ∑numRands = 0
    ∑weighted = 0
    
    for i in (2 : length(numTets_by_index))
        ∑numRands +=   numTets_by_index[i]
        ∑weighted += i*numTets_by_index[i]
    end


    ∑weighted / 
    ∑numRands
end

# ╔═╡ e9cd89ce-f333-11ea-3336-7f5d3574d63d
numTets = MCe_∑rands!([], 100)

# ╔═╡ fab8359a-f333-11ea-3147-958990062899
e_est_from_tets(numTets)

# ╔═╡ 0552c218-f334-11ea-3472-b71f496c5234
MCe_∑rands!(numTets, 900)

# ╔═╡ 624a62c8-f334-11ea-2847-6d53423efd65
e_est_from_tets(numTets)

# ╔═╡ 72f6f8f2-f334-11ea-1858-f3fccd3a3c0b
MCe_∑rands!(numTets, 9000)

# ╔═╡ 7c13f9da-f334-11ea-2ad6-59453015db59
e_est_from_tets(numTets)

# ╔═╡ 0efddefe-f336-11ea-3eb3-e5e6d92551a7
MCe_∑rands!(numTets, 90_000)

# ╔═╡ 18366024-f336-11ea-2a3b-0d8f36ba850e
e_est_from_tets(numTets)

# ╔═╡ ae9cd56c-f336-11ea-2b96-7307425123d3
MCe_∑rands!(numTets, 900_000)

# ╔═╡ b705b9d2-f336-11ea-1daa-0f488ca41936
e_est_from_tets(numTets)

# ╔═╡ 7635f79e-f338-11ea-0072-8936cce52f04
MCe_∑rands!(numTets, 9_000_000)

# ╔═╡ 8a3e2c5c-f338-11ea-3f4c-95b818a17841
e_est_from_tets(numTets)

# ╔═╡ 87243e50-f70f-11ea-1066-eb6f4c4967a1
function MCe_darts(rows::Int, cols::Int, width, height) 
	
    hits_by_row_col = zeros(Int8, rows, cols) 
	hits_coords = []
   
    for dart in (1 : rows*cols)
		rnd = rand()
		y = height*rnd
		row = convert(Int64, ceil(rows*rnd))
		
		
		rnd = rand()
		x = width*rnd
		col = convert(Int64, ceil(cols*rnd))	
		
		hits_by_row_col[row, col] += 1
		
		push!(hits_coords, (x, y))		
    end
    
	
    hits_by_row_col, hits_coords
end

# ╔═╡ 869e3ffa-f712-11ea-0e66-c30e457600d5
hits_rowcol, hits_coords = MCe_darts(10, 10, 100, 100)

# ╔═╡ 88a5095c-f713-11ea-2fff-3f7b8b12ec6d
hit_cnts2 = hit_cnts_from_targets(hits_rowcol)

# ╔═╡ Cell order:
# ╠═ec55d1a4-f2a4-11ea-21ba-d5d9e9fab3e5
# ╠═30ef477e-f33a-11ea-257f-fdb9ba8bf110
# ╠═fa15ffc8-f2a2-11ea-39cd-27af9b2559d0
# ╠═171ee65e-f33c-11ea-012e-3b4e1ba23fc9
# ╠═97cf78d8-f33a-11ea-17fa-6157a0fad0b1
# ╠═3eaca8a0-f33c-11ea-14f7-258438f5b563
# ╟─750ae004-f339-11ea-2fdc-3386c26e35a0
# ╠═a7ea64ae-f33e-11ea-00ba-d7e2dd0003ce
# ╟─b38f332c-f328-11ea-3ab7-c5ed3b530bbb
# ╟─d5a282a6-f333-11ea-209f-2d6d48671f3e
# ╠═e9cd89ce-f333-11ea-3336-7f5d3574d63d
# ╠═fab8359a-f333-11ea-3147-958990062899
# ╠═0552c218-f334-11ea-3472-b71f496c5234
# ╠═624a62c8-f334-11ea-2847-6d53423efd65
# ╠═72f6f8f2-f334-11ea-1858-f3fccd3a3c0b
# ╠═7c13f9da-f334-11ea-2ad6-59453015db59
# ╠═0efddefe-f336-11ea-3eb3-e5e6d92551a7
# ╠═18366024-f336-11ea-2a3b-0d8f36ba850e
# ╠═ae9cd56c-f336-11ea-2b96-7307425123d3
# ╠═b705b9d2-f336-11ea-1daa-0f488ca41936
# ╠═7635f79e-f338-11ea-0072-8936cce52f04
# ╠═8a3e2c5c-f338-11ea-3f4c-95b818a17841
# ╠═87243e50-f70f-11ea-1066-eb6f4c4967a1
# ╠═869e3ffa-f712-11ea-0e66-c30e457600d5
# ╠═88a5095c-f713-11ea-2fff-3f7b8b12ec6d
