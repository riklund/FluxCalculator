#!/bin/bash
function Main()
{
	InitRLlib
	if [ -f "makefile" ]
	then
		chmod a+w makefile
	fi
	cp makefile.in makefile
	#To prevent accidental modification of makefile instead of makefile.in
	chmod a-w makefile 
	echo "Setup completed without errors."
}

#Run setup script on RLlib
function InitRLlib()
{
	if [ ! -f "RLlib/setup.sh" ]
	then
		echo "Could not find RLlib/setup.sh, fatal error."
		exit 4
	fi
	cd RLlib
	echo -n "RLlib : "
	./setup.sh
	if [ "$?" -ne "0" ]
	then
		echo "RLlib setup returned nonzero exit code, fatal error."
		exit 5
	fi
	cd ..
}






#Run the main function.
Main "$*"