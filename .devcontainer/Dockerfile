FROM ubuntu:oracular

# Update system
RUN apt-get update -y && apt-get upgrade -y

# Install general dev packages
RUN apt-get install -y \
    apt-utils cmake g++ clang-format git curl wget sudo

# Install Sparselizard dependencies
RUN apt-get install -y \
    libopenblas-dev \
    libopenmpi-dev \
    libmumps-dev \
    libmetis-dev \
    petsc-dev \
    slepc-dev \
    libgmsh-dev \
    libomp-dev

# Use the provided ubuntu user instead of root
ARG USERNAME=ubuntu
RUN echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME
USER $USERNAME

# Install Oh My Zsh for convenience
RUN sh -c "$(wget -O- https://github.com/deluan/zsh-in-docker/releases/download/v1.2.1/zsh-in-docker.sh)" -- \
    -t robbyrussell \
    -p git

# Set the default shell to zsh
RUN sudo chsh -s $(which zsh)
CMD ["zsh"]