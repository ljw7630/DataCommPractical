set num=%1


set srcx=%2


set "str="import receiver_preprocess; receiver_preprocess.getRandomPairs(%num%, %srcx%)""

python -c %str%
