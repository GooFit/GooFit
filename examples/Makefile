.PHONY: all clean

all: 2d_plot/2d_plot addition/addition chisquare/chisquare convolution/convolution dalitz/dalitz DP4/DP4 exponential/exponential pipipi0DPFit/pipipi0DPFit product/product SigGen/SigGen simpleFit/simpleFit TDDP4/TDDP4 zachFit/zachFit

2d_plot/2d_plot: 2d_plot/2d_plot.cu
	cd 2d_plot && $(MAKE)

addition/addition: addition/addition.cu
	cd addition && $(MAKE)	

chisquare/chisquare: chisquare/chisquare.cu
	cd chisquare && $(MAKE)

convolution/convolution: convolution/convolution.cu
	cd convolution && $(MAKE)

dalitz/dalitz: dalitz/dalitz.cu
	cd dalitz && $(MAKE)

DP4/DP4: DP4/DP4.cu
	cd DP4 && $(MAKE)

exponential/exponential: exponential/exponential.cu
	cd exponential && $(MAKE)

pipipi0DPFit/pipipi0DPFit: pipipi0DPFit/pipipi0DPFit.cu
	cd pipipi0DPFit && $(MAKE)

product/product: product/product.cu
	cd product && $(MAKE)

SigGen/SigGen: SigGen/SigGen.cu
	cd SigGen && $(MAKE)

simpleFit/simpleFit: simpleFit/simpleFit.cu
	cd simpleFit && $(MAKE)

TDDP4/TDDP4: TDDP4/TDDP4.cu
	cd TDDP4 && $(MAKE)

zachFit/zachFit: zachFit/zachFit.cu
	cd zachFit && $(MAKE)

clean:
	cd 2d_plot && $(MAKE) clean
	cd addition && $(MAKE) clean
	cd chisquare && $(MAKE) clean
	cd convolution && $(MAKE) clean
	cd dalitz && $(MAKE) clean
	cd DP4 && $(MAKE) clean
	cd exponential && $(MAKE) clean
	cd pipipi0DPFit && $(MAKE) clean
	cd product && $(MAKE) clean
	cd SigGen && $(MAKE) clean
	cd simpleFit && $(MAKE) clean
	cd TDDP4 && $(MAKE) clean
	cd zachFit && $(MAKE) clean
