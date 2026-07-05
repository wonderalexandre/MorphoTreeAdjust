#pragma once

#include "../details/NonIncreasingCasfKernel.hpp"

namespace debug::non_increasing {

struct NonIncreasingCasfDiagnosticsOptions {
    bool collectMergeEvents = false;
    bool collectPhaseOracleComparisons = false;
    bool phaseCandidateEquivalenceGuard = false;
    bool stopOnFirstInvariantMismatch = false;
    bool stopOnFirstPhaseRebuildFallback = false;
    bool stopOnFirstPhaseReleaseRejection = false;
    bool collectFirstInvariantMismatchOracle = false;
    bool leafExplicitFastPath = true;
    bool mcpcDirectPreflightDiagnostics = false;
    bool mcpcDirectV1StrictPrecheck = false;
    bool mcpcNoCopyDirectCommit = false;
    bool mcpcValidatedDirectCommit = false;
    bool operationalMergeEventInstrumentation = true;
    int maxAppliedMergeEventsPerThreshold = 0;
    long long maxThresholdMicros = 0;
    bool captureLocalOracle = false;
    float localOracleThreshold = 0.0f;
    ::morphotree::details::non_increasing::NonIncreasingPhase localOraclePhase =
        ::morphotree::details::non_increasing::NonIncreasingPhase::MaxToMin;
    NodeId localOracleNodeId = InvalidNode;
    ::morphotree::details::non_increasing::NonIncreasingStrippedLocalUpdateOptions strippedLocalUpdateOptions;
    ::morphotree::details::non_increasing::NonIncreasingStrippedCommitMode strippedCommitMode =
        ::morphotree::details::non_increasing::NonIncreasingStrippedCommitMode::SafeCopy;
};

inline void applyNonIncreasingCasfDiagnostics(
    ::morphotree::details::non_increasing::NonIncreasingCasfKernel &runner,
    const NonIncreasingCasfDiagnosticsOptions &options) {
    runner.setStrippedLocalUpdateOptions(options.strippedLocalUpdateOptions);
    runner.setStrippedCommitMode(options.strippedCommitMode);
    runner.setLeafExplicitFastPathEnabled(options.leafExplicitFastPath);
    runner.setPhaseCandidateEquivalenceGuardEnabled(options.phaseCandidateEquivalenceGuard);
    runner.setMcpcDirectPreflightDiagnosticsEnabled(
        options.mcpcDirectPreflightDiagnostics);
    runner.setMcpcDirectV1StrictPrecheckEnabled(options.mcpcDirectV1StrictPrecheck);
    runner.setMcpcNoCopyDirectCommitEnabled(options.mcpcNoCopyDirectCommit);
    runner.setMcpcValidatedDirectCommitEnabled(options.mcpcValidatedDirectCommit);
    runner.setOperationalMergeEventInstrumentationEnabled(
        options.operationalMergeEventInstrumentation);
    runner.setMergeEventCollectionEnabled(options.collectMergeEvents);
    runner.setPhaseOracleComparisonEnabled(options.collectPhaseOracleComparisons);
    runner.setStopOnFirstInvariantMismatchEnabled(options.stopOnFirstInvariantMismatch);
    runner.setStopOnFirstPhaseRebuildFallbackEnabled(
        options.stopOnFirstPhaseRebuildFallback);
    runner.setStopOnFirstPhaseReleaseRejectionEnabled(
        options.stopOnFirstPhaseReleaseRejection);
    runner.setFirstInvariantMismatchOracleCollectionEnabled(
        options.collectFirstInvariantMismatchOracle);
    runner.setMaxAppliedMergeEventsPerThreshold(
        options.maxAppliedMergeEventsPerThreshold);
    runner.setMaxThresholdMicros(options.maxThresholdMicros);
    if (options.captureLocalOracle) {
        runner.setLocalOracleCaptureTarget(
            options.localOracleThreshold,
            options.localOraclePhase,
            options.localOracleNodeId);
    } else {
        runner.clearLocalOracleCaptureTarget();
    }
}

} // namespace debug::non_increasing
