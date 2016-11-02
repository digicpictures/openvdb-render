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

ProgressBar::ProgressBar(const MString& msg, const uint32_t max_progress, const bool is_interruptable)
  : m_show_progress_bar(MGlobal::mayaState() == MGlobal::kInteractive),
    m_is_interruptable(is_interruptable)
{
    if (!m_show_progress_bar)
        return;

    reset(msg, max_progress);
}

ProgressBar::~ProgressBar()
{
    if (!m_show_progress_bar)
        return;

    endProgress();
}

void ProgressBar::reset(const MString& msg, const uint32_t max_progress)
{
    if (!m_show_progress_bar)
        return;

    m_max_progress = max_progress;
    m_progress = 0;
    beginProgress(msg);
}

void ProgressBar::addProgress(uint32_t progress_to_add)
{
    if (!m_show_progress_bar)
        return;

    const auto old_progress = m_progress.fetch_and_add(progress_to_add);
    const auto old_percents = uint8_t(old_progress * 100 / m_max_progress);
    const auto new_percents = uint8_t((old_progress + progress_to_add) * 100 / m_max_progress);
    if (new_percents > old_percents) {
        tbb::mutex::scoped_lock lock(m_mutex);
        MGlobal::executeCommand(format("progressBar -e -pr ^1s $gMainProgressBar", new_percents));
    }
}

void ProgressBar::setProgress(const int percent)
{
    if (!m_show_progress_bar)
        return;

    MGlobal::executeCommand(format("progressBar -e -pr ^1s $gMainProgressBar; refresh", percent));
}

bool ProgressBar::isCancelled() const
{
    if (!m_show_progress_bar || !m_is_interruptable)
        return false;

    tbb::mutex::scoped_lock lock(m_mutex);

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
    MGlobal::executeCommand(format("progressBar -e -bp -ii ^1s -st \"^2s\" -max 100 $gMainProgressBar", m_is_interruptable, msg));
}

void ProgressBar::endProgress() const
{
    MGlobal::executeCommand("progressBar -e -ep $gMainProgressBar");
}
