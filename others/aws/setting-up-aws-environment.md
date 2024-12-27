# setting up aws environment

## Setting up docker

```bash
# setting up aws environment <https://docs.docker.com/engine/install/ubuntu/>
# sudo issues <https://www.baeldung.com/linux/docker-permission-denied-daemon-socket-error>
# Add Docker's official GPG key:
sudo apt-get update
sudo apt-get install ca-certificates curl gnupg
sudo install -m 0755 -d /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
sudo chmod a+r /etc/apt/keyrings/docker.gpg

# Add the repository to Apt sources:
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$VERSION_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

# setting up docker
sudo chmod 666 /var/run/docker.sock
sudo docker run -it ubuntu:latest /bin/bash
sudo groupadd docker
getent group docker
awk -F':' '/docker/{print $4}' /etc/group
sudo usermod -aG docker ubuntu
getent group docker
awk -F':' '/docker/{print $4}' /etc/group
sudo systemctl restart docker.service
sudo systemctl status docker.service
docker run -it ubuntu:latest /bin/bash # check

# setting up java <https://www.theserverside.com/blog/Coffee-Talk-Java-News-Stories-and-Opinions/How-do-I-install-Java-on-Ubuntu>
# <https://www.rosehosting.com/blog/how-to-install-java-se-19-on-ubuntu-22-04/>
sudo apt-get update -y && sudo apt-get upgrade -y
sudo apt install openjdk-19-jdk openjdk-19-jre

# setting path
export PATH="/home/ubuntu/.local/bin:$PATH"
```

