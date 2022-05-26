# EchoProber-server

MATLAB server for real-time audio visualization and analysis.

Check out https://github.com/felixnie/EchoProber for the client app for Android.

**Please upgrade to MATLAB R2021b or later to use tcpserver.**

## Functions

1. Manage connections from multiple clients.
2. Continuous data receiving and plotting.
3. Manage multiple figures: figures will refresh in background instead of popping up.
4. Resolve message from clients: short message as commands, long message as recorded data.
5. Remote control on clients: short message as commands, long message as chirp data.

## To-do

1. (Pending) Parallel processing for analysis tasks.
2. (Pending) Deal with the lagging when the collected data array is large.

## Notes

This server runs in the background. 
One can run other scripts to communicate with clients once the server is set up. 
The MATLAB live scripts are some examples for remote control and data collection.

Before restarting the server, please disconnect all the clients and check the host IP and ports.
