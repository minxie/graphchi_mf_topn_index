export GRAPHCHI_ROOT=../
for m in '50'
do
	for l in  'le-3'
	do
		../bin/myapps/sgd --D=$m --training=../../data/yahoomusic/ratings_kddcup_mm --tol= 1e-1 --result=../../result/gc_yahoomusic_res \
		  --sgd_gamma=1e-4 --sgd_lambda=1e-4 --max_iter=100 --n_top=50 --quiet=1  --ff=4 --tol=1e-1
	done
done
