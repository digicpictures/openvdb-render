#pragma once

#include <mutex>

#include <maya/MGlobal.h>
#include <maya/MString.h>

// Based on the gpuCache plugin from the devkit.

class ProgressBar
{
public:
    ProgressBar(const MString& msg, const bool is_interruptable = false);
    ~ProgressBar();

    void reset(const MString& msg);

    // The public methods below are thread-safe.
    void stepOnePercent() const;
    void setProgress(const int percent);
    bool isCancelled() const;

private:
    ProgressBar(const ProgressBar&) = delete;
    ProgressBar& operator=(const ProgressBar&) = delete;
    ProgressBar(ProgressBar&&) = delete;
    ProgressBar&& operator=(ProgressBar&&) = delete;

    void beginProgress(const MString& msg) const;
    void endProgress() const;

    bool m_show_progress_bar;
    bool m_is_interruptable;
    mutable std::mutex m_mutex;
};
