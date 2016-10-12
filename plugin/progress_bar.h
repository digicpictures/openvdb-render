#pragma once

#include <maya/MGlobal.h>
#include <maya/MString.h>

class ProgressBar
{
public:
    ProgressBar(const MString& msg, unsigned int max)
    {
        // Display a progress bar if Maya is running in UI mode
        fShowProgress = (MGlobal::mayaState() == MGlobal::kInteractive);
        reset(msg, max);
    }
    void reset(const MString& msg, unsigned int max)
    {
        MStatus status;
        beginProgress(msg, max);
    }
    ~ProgressBar()
    {
        endProgress();
    }
    void stepProgress() const
    {
        if (fShowProgress) {
            MGlobal::executeCommand("progressBar -e -s 1 $gMainProgressBar");
        }
    }
    bool isCancelled() const
    {
        int isCancelled = 0;
        if (fShowProgress) {
            MGlobal::executeCommand("progressBar -q -ic $gMainProgressBar", isCancelled);
        }
        if (isCancelled) {
            MStatus status;
            const MString interruptMsg = "Interrupted by user";
            MGlobal::displayInfo(interruptMsg);
            return true;
        }
        return false;
    }
private:
    // Forbidden and not implemented.
    ProgressBar(const ProgressBar&);
    const ProgressBar& operator=(const ProgressBar&);
    void beginProgress(const MString& msg, unsigned int max) const
    {
        if (fShowProgress) {
            MString maxValue, progressBarCmd;
            // Progress from 0 to max
            if (max <= 0) {
                max = 1;
            }
            maxValue += max;
            // Clear previous isCancelled flag
            MGlobal::executeCommand("progressBar -e -bp -ii 1 $gMainProgressBar");
            MGlobal::executeCommand("progressBar -e -ep $gMainProgressBar");
            // Initialize the progress bar
            progressBarCmd.format("progressBar -e -bp -ii 1 -st \"^1s\" -max ^2s $gMainProgressBar",
                msg, maxValue);
            MGlobal::executeCommand(progressBarCmd);
        }
    }
    void endProgress() const
    {
        if (fShowProgress) {
            MGlobal::executeCommand("progressBar -e -ep $gMainProgressBar");
        }
    }
    bool fShowProgress;  // whether to show the progress bar
};

