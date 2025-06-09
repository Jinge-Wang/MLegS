# Connect to island.me.berkeley.edu
# Public key: key_Rsync for rsync use only
rsync -auvz -e 'ssh -p 7777 -i /home/jinge/.ssh/key_Rsync' /home/jinge/Desktop/InitialValueCode jinge@island.me.berkeley.edu:/home/jinge/storage
