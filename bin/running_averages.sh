#! /bin/bash

if [ "$#" -ge "1" ]; then
	if [ -e $1 ]; then
		if [ $# -eq "1" ]; then
			BIN="10"
		else
			BIN=$2
		fi
		awk -v bin=$BIN '{ 
			num++; 
			x[num]=$1; 
			y[num]=$2; 
			ye[num]=$3; 
		}END{ 
			if (bin>0) {
				for ( i=1;i<=bin;i++ ) {
					xsum+=x[i]; 
					ysum+=y[i]; 
					y2sum+=y[i]^2; 
				} 
				xav = xsum/bin ;
				yav = ysum/bin ;
				ystd =  sqrt(bin/(bin-1)) * sqrt( y2sum/bin - (ysum/bin)^2 ) ;
				yerr =  ystd / sqrt(bin) ;
				printf( "%f   %.17e   %.17e   %.17e\n", xav, yav, yerr, ystd );
				for(j=1;j<=num-bin;j++) {
					xsum=xsum+x[bin+j]-x[j]; 
					ysum=ysum+y[bin+j]-y[j]; 
					y2sum=y2sum+y[bin+j]^2-y[j]^2; 
					xav = xsum/bin ;
					yav = ysum/bin ;
					ystd =  sqrt(bin/(bin-1)) * sqrt( y2sum/bin - (ysum/bin)^2 ) ;
					yerr =  ystd / sqrt(bin) ;
					printf( "%f   %.17e   %.17e   %.17e\n", xav, yav, yerr, ystd );
				} 
			}
			else {
				for ( i=1;i<=num;i++ ) {
					xsum+=x[i]; 
					ysum+=y[i]; 
					y2sum+=y[i]^2; 
				} 
				xav = xsum/num ;
				yav = ysum/num ;
				ystd =  sqrt(num/(num-1)) * sqrt( y2sum/num - (ysum/num)^2 ) ;
				yerr =  ystd / sqrt(num) ;
				#printf( "%f   %.17e   %.17e   %.17e\n", xav, yav, yerr, ystd );
				printf( "0.   %.17e   %.17e   %.17e\n", yav, yerr, ystd );
			}
			}' $1
	else
		echo "ERROR: file $1 do not exist!!!"
	fi
else
	echo "WRITE: running_averages.sh <file:x_y_dy> <bin>"
fi


