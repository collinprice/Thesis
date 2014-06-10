sudo apt-get -y update
sudo apt-get -y install vim ifeffit make gcc g++

if [ ! -f /vagrant/EXAFS-OEC/output.dcd.full ]; then
    cd /vagrant/EXAFS-OEC 
	wget https://dl.dropboxusercontent.com/u/639928/thesis/output.dcd.full
fi

if ! which namd2 >/dev/null; then
	cd /tmp
	wget https://dl.dropboxusercontent.com/u/639928/thesis/NAMD_2.9_Linux-x86.tar.gz
	tar -zxvf NAMD_2.9_Linux-x86.tar.gz
	rm NAMD_2.9_Linux-x86.tar.gz
	cd NAMD_2.9_Linux-x86
	sudo cp -r * /usr/local/bin/
	cd ..
	rm -rf NAMD_2.9_Linux-x86
fi

if ! which vmd >/dev/null; then
	cd /tmp
	wget https://dl.dropboxusercontent.com/u/639928/thesis/vmd-1.9.1.bin.LINUX.opengl.tar.gz
	tar -zxvf vmd-1.9.1.bin.LINUX.opengl.tar.gz vmd
fi

if ! which vmd >/dev/null; then

	sudo apt-get -y install libglu1-mesa libxinerama-dev libxft-dev libxi-dev

	cd /tmp
	wget https://dl.dropboxusercontent.com/u/639928/thesis/vmd-1.9.1.bin.LINUX.opengl.tar.gz
	tar -zxvf vmd-1.9.1.bin.LINUX.opengl.tar.gz
	rm vmd-1.9.1.bin.LINUX.opengl.tar.gz
	cd vmd-1.9.1

	./configure LINUX
	cd src
	sudo make install

	cd ../../
	rm -rf vmd-1.9.1
fi