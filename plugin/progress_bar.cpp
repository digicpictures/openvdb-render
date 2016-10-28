#include "progress_bar.h"

namespace {

    template <typename T>
    inline MString toMString(const T& value)
    {
        MString res;
        res.set(value);
        return res;
    }

    template <>
    inline MString toMString<MString>(const MString& value)
    {
        return value;
    }

    template <typename... Args>
    inline MString format(const MString& format, Args&&... args)
    {
        MString res;
        res.format(format, toMString(std::forward<Args>(args))...);
        return res;
    }

} // unnamed namespace

ProgressBar::ProgressBar(const MString& msg, const bool is_interruptable)
  : m_show_progress_bar(MGlobal::mayaState() == MGlobal::kInteractive),
    m_is_interruptable(is_interruptable)
{
    if (!m_show_progress_bar)
        return;

    reset(msg);
}

ProgressBar::~ProgressBar()
{
    if (!m_show_progress_bar)
        return;

    endProgress();
}

void ProgressBar::reset(const MString& msg)
{
    if (!m_show_progress_bar)
        return;

    beginProgress(msg);
}

void ProgressBar::stepOnePercent() const
{
    if (!m_show_progress_bar)
        return;

    std::lock_guard<std::mutex> lock(m_mutex);
    MGlobal::executeCommand("progressBar -e -s 1 $gMainProgressBar; refresh");
}

void ProgressBar::setProgress(const int percent)
{
    if (!m_show_progress_bar)
        return;

    std::lock_guard<std::mutex> lock(m_mutex);

    //MString cmd, percent_str;
    //percent_str.set(percent);
    //cmd.format("progressBar -e -pr ^1s $gMainProgressBar", percent_str);
    std::cout << "ProgressBar::setProgress(" << percent << ")" << std::endl;
    MGlobal::executeCommand(format("progressBar -e -pr ^1s $gMainProgressBar; refresh", percent));
}

bool ProgressBar::isCancelled() const
{
    if (!m_show_progress_bar || !m_is_interruptable)
        return false;

    std::lock_guard<std::mutex> lock(m_mutex);

    int is_cancelled = 0;
    MGlobal::executeCommand("progressBar -q -ic $gMainProgressBar", is_cancelled);
    if (is_cancelled)
        MGlobal::displayInfo("Interrupted by user");

    return is_cancelled != 0;
}

void ProgressBar::beginProgress(const MString& msg) const
{
    // Clear previous isCancelled flag.
    MGlobal::executeCommand("progressBar -e -bp -ii 1 $gMainProgressBar");
    MGlobal::executeCommand("progressBar -e -ep $gMainProgressBar");
    // Initialize the progress bar.
    MString cmd;
    cmd.format("progressBar -e -bp -ii 1 -st \"^1s\" -max 100 $gMainProgressBar", msg);
    MGlobal::executeCommand(format("progressBar -e -bp -ii ^1s -st \"^2s\" -max 100 $gMainProgressBar", m_is_interruptable, msg));
}

void ProgressBar::endProgress() const
{
    MGlobal::executeCommand("progressBar -e -ep $gMainProgressBar");
}
