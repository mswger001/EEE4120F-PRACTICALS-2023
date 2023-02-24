
#Author: ?????

using Pkg

#_____________Uncomments the packages below once they have been installed________
using WAV
using Plots
using TickTock
using Statistics
using CSV
using DataFrames
#using StatsPlots

#this function installs all packages needed to run this script
function installFunc()
	Pkg.add("WAV")
	Pkg.add("StatsPlots")
	Pkg.add("Plots")
	Pkg.add("TickTock")
	Pkg.add("Statistics")
	Pkg.add("DataFrames")
    Pkg.add("CSV")
	print("------Packages Installed---------")
	print(Pkg.status())
end


function main()
	#installFunc() #comment this out once you have run it once and all packages have been installed
	whiteNoise = (rand(48000)*2).-1

	WAV.wavwrite(whiteNoise, "whiteNoise.wav", Fs=4800) #sample freq is 4800Hz
	
	tick()
	whiten = createwhiten(10); #this will create 1000 seconds of white noise
	print(tok())
	
	WAV.wavwrite(whiten , "white_noise_sound2.wav", Fs=4800)
	i =  Plots.histogram(whiteNoise, title= "histogram plot for white noise", xlabel ="number generated",ylabel="frequency")
	
	h = Plots.histogram(whiten, title= "histogram plot for whiten noise", xlabel ="number generated",ylabel="frequency")
	
	Plots.display(h)
	savefig(h,"whitenHistogram.png")
	savefig(i,"whiteNoiseHist.png")
	readline() #this will stop the program at this point till you press enter

	createwhitenExperiment()
	correlation_experiment()
	corrShiftedSignals()

	
end

function correlation_experiment()
	
	sample_siz = []
	custom_corr_valu =  []
	custom_corr_tim = []
	statistics_cor_valu = []
	statistics_cor_tim  = []

	custom_corr_value =  0
	custom_corr_time = 0
	statistics_cor_value = 0
	statistics_cor_time  =0

	for j in 1:20
		sample_size = 100 
		for i in 1:10
			whiten = createwhiten(sample_size); #this will create 1000 seconds of white noise
			tick()
			custom_corr_value = corr(whiten,whiten)
			custom_corr_time = tok()
				
		
			tick()
			statistics_cor_value = Statistics.cor(whiten,whiten)
			statistics_cor_time = tok()
		
			push!(sample_siz, sample_size)
			push!(custom_corr_valu, custom_corr_value)
			push!(custom_corr_tim, custom_corr_time)
			push!(statistics_cor_valu, statistics_cor_value)
			push!(statistics_cor_tim, statistics_cor_time)
			
			

			sample_size = sample_size + 100;
		end

	end


	df = DataFrame(
		sample_size = sample_siz,
		custom_corr_value = custom_corr_valu,
		custom_corr_time = custom_corr_tim,
		statistics_corr_value = statistics_cor_valu,
		statistics_corr_time = statistics_cor_tim,
	)
	dir = joinpath(pwd(), "correlationExperiment.csv")
	CSV.write(dir, df)
	
end

function createwhitenExperiment()
	sample_size = 0
	sz =[]
	whiteNoiseTT =[]
	whitenNoiseTT =[]
	for j in 1:20
		sample_size = 100 
		for i in 1:10
			tick() # we taking the time it starts
			whiteNoise = (rand(4800*sample_size)*2).-1 
			time_taken = tok()

			tick() # we taking the time it starts
			whiten = createwhiten(sample_size); #this will create 10 seconds of white noise
			whiten_time_taken = tok()

			push!(sz,sample_size)
			push!(whitenNoiseTT,time_taken)
			push!(whiteNoiseTT,whiten_time_taken)	

			sample_size = sample_size+100

		end
	end
	df = DataFrame(
		sample_size = sz,
		time_taken_white_noise = whiteNoiseTT,
		time_taken_whiten_noise = whitenNoiseTT,
 	)
	dir = joinpath(pwd(), "noise_creation_experiment.csv")
	CSV.write(dir, df)

end


function createwhiten(N)
	#f = 1/t where in this case we want
	# freq to always be 48khz
	# if for 1 sec  we use 48000 freq then we sample at the same rate for N times
	arr = Array{Float64,1}()
	j=0
	n=1
	max =N*48*100
	for i in 1:max
		half = max/2
		if i>half 
			n= (-1)
		end	
		push!(arr,  float(rand())*n)
		j=i
	end
	#print("Number of samples: (")
	#print(j)
	#print(",)")
	return arr
end 

function corr(x,y)
	sumX = sum(x)
	sumY = sum(y)
	avgX = sumX/length(x)
	avgY = sumY/length(y)

	numerator = 0
	denominator = 0
	diffx2summ = 0
	diffy2summ = 0
	for i in 1:length(x)
		difx = x[i]- avgX
		dify = y[i] - avgY
		mulXY = difx*dify
		numerator = numerator+mulXY
		difx2 = difx*difx
		diffx2summ = diffx2summ +difx2
		dify2 = dify*dify
		diffy2summ = diffy2summ +dify2
	end
	#denominator = (1/(diffx2summ*diffx2summ))*(1/(diffy2summ*diffy2summ))
	
	result = numerator/sqrt(diffx2summ*diffy2summ)
	return result
end 

function corrShiftedSignals()
	
	sinewaves = []
	sz = []
	fr =[]
	shift = []
	count = 0
	sample_size = 100 
	for j in 1:3
		for k in 1:3			
			freq = k*5	
			for inn in 0:3	
				count = count+1			
				t = range(0, stop=1, length=sample_size)
				sine_wave = sin.(2π * freq .* t)
				shifted_sine_wave = circshift(sine_wave,inn*5)
				push!(shift, (inn*5))
				push!(sinewaves,shifted_sine_wave)
				push!(sz,sample_size)	
				push!(fr,freq)		
			end
		end
		sample_size = sample_size*10
	end

	#sample_size,freq,shift,correlation with the 0 shift waveof same sample size,
	corr_value =[]
	unshiftSinewave1000 =0
	unshiftSinewave100 =0
	unshiftSinewave10000 =0
	for i in 1:count
		outname= string(fr[i])*"_"*string(sz[i])*"_"*string(shift[i])*".png"
		plt = Plots.plot(sinewaves[i], xlabel="time", ylabel="Y")
		savefig(plt,outname)
		#only have corr for the same sample_siz
		if sz[i]==100			
			if  shift[i]==0
				unshiftSinewave100 = sinewaves[i] 
			end
			if unshiftSinewave100==0
				continue
			else	
				push!(corr_value,Statistics.cor(unshiftSinewave100,sinewaves[i]))

			end			
		end
		if sz[i]==1000		
			if  shift[i]==0
				unshiftSinewave1000 = sinewaves[i] 
			end
			if unshiftSinewave1000==0
				continue
			else
					
				push!(corr_value,Statistics.cor(unshiftSinewave1000,sinewaves[i]))

			end			
		end
		if sz[i]==10000
			if  shift[i]==0
				unshiftSinewave10000 = sinewaves[i] 
			end
			if unshiftSinewave10000==0
				continue
			else	
				push!(corr_value,Statistics.cor(unshiftSinewave10000,sinewaves[i]))

			end			
		end
	end
	df = DataFrame(
		sample_size = sz,
		frequency = fr,
		shift  = shift,
		corrToZero = corr_value,)

	dir = joinpath(pwd(), "shifted_sine_wave_experiment.csv")
	CSV.write(dir, df)

	
	
		




end
	

main()
