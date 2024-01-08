#!/bin/bash
if [[ "$OSTYPE" == "linux-gnu" ]]; then
  shell_profile=$HOME/.bashrc
  os_opts=''
elif [[ "$OSTYPE" == "darwin"* ]]; then
  shell_profile=$HOME/.bash_profile
  os_opts='-v/Users:/Users:delegated'
else
  shell_profile=$HOME/.bashrc
  os_opts=''
fi

echo ""
echo "Hello!.. this installs 3 lines to your .bashrc (thingy that is run when you start a Terminal)"
echo "These 3 lines setup an 'osiris' command that uses Docker to do compile+run Osiris!"
echo ""
echo "     Adding the 3 lines...."
echo 'osiris() {' >> $shell_profile
echo '    docker run -it -v$(pwd):/osiris:delegated' "$os_opts" '--user $(id -u):$(id -g) -w="/osiris" atableman/osirispic $1 $2 $3 $4 $5 $6' >> $shell_profile
echo '}' >> $shell_profile
echo ""
echo "     Done.... Restart (ie restart your Terminal/log out then back into Terminal to see the effects)"
echo ""


