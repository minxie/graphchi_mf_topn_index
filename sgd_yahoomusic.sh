export GRAPHCHI_ROOT=../
for m in '10' '20' '30' '40' '50'
do
	for l in  'le-3'
	do
		../bin/myapps/sgd2 --D=$m --training=../../../../../../../research/connections/people/minxie/Data/KDDCup11/track1/ratings_kddcup_alt_mm --result=../../result/gc_yahoomusic_res \
		  --sgd_gamma=1e-6 --sgd_lambda=1e-6 --max_iter=100 --n_top=50 --quiet=1  --tol=1e-1
	done
done
