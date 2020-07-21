getenv			= true
universe		= vanilla
+RequiresWholeMachine	= True
Requirements		= isapollo == True
notification    	= Always
notify_user     	= hocamachoc@gmail.com
executable     		= script/driver.sh
output         		= log/test.$(ClusterId).$(ProcId).out
error          		= log/test.$(ClusterId).$(ProcId).err
log            		= log/test.$(ClusterId).$(ProcId).log
queue
