package edu.upenn.seas.precise.verisig.utils;

import java.io.*;

/**
 * This class is intended to be used with the SystemCommandExecutor
 * class to let users execute system commands from Java applications.
 *
 * This class is based on work that was shared in a JavaWorld article
 * named "When System.exec() won't". That article is available at this
 * url:
 *
 * http://www.javaworld.com/javaworld/jw-12-2000/jw-1229-traps.html
 *
 * Documentation for this class is available at this URL:
 *
 * http://devdaily.com/java/java-processbuilder-process-system-exec
 *
 *
 * Copyright 2010 alvin j. alexander, devdaily.com.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser Public License for more details.

 * You should have received a copy of the GNU Lesser Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Please ee the following page for the LGPL license:
 * http://www.gnu.org/licenses/lgpl.txt
 *
 */
class ThreadedStreamHandler extends Thread
{
    private InputStream inputStream;
    private StringBuilder outputBuffer = new StringBuilder();

    /**
     * A simple constructor for when the sudo command is not necessary.
     * This constructor will just run the command you provide, without
     * running sudo before the command, and without expecting a password.
     *
     * @param inputStream
     */
    ThreadedStreamHandler(InputStream inputStream)
    {
        this.inputStream = inputStream;
    }

    public void run()
    {
        try (BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(inputStream))) {
            String line = null;
            while ((line = bufferedReader.readLine()) != null) {
                outputBuffer.append(line + "\n");
            }
        } catch (Throwable ioe) {
            // TODO handle this better
            ioe.printStackTrace();
        }
        // ignore this one
    }

    public StringBuilder getOutputBuffer()
    {
        return outputBuffer;
    }

}








