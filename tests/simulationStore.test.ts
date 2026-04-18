import { act } from 'react'
import { useSimulationStore } from '../src/state/simulationStore'

describe('simulation store', () => {
  beforeEach(() => {
    useSimulationStore.setState(useSimulationStore.getState(), true)
    useSimulationStore.getState().reset()
    useSimulationStore.setState({ simSpeed: 1, isPlaying: false })
  })

  it('toggles play state', () => {
    expect(useSimulationStore.getState().isPlaying).toBe(false)
    act(() => useSimulationStore.getState().play())
    expect(useSimulationStore.getState().isPlaying).toBe(true)
    act(() => useSimulationStore.getState().pause())
    expect(useSimulationStore.getState().isPlaying).toBe(false)
  })

  it('updates time and speed', () => {
    act(() => useSimulationStore.getState().setTime(10))
    act(() => useSimulationStore.getState().setSimSpeed(3))
    expect(useSimulationStore.getState().time).toBe(10)
    expect(useSimulationStore.getState().simSpeed).toBe(3)
  })

  it('toggles overlays', () => {
    const initial = useSimulationStore.getState().showWind
    act(() => useSimulationStore.getState().toggleOverlay('showWind'))
    expect(useSimulationStore.getState().showWind).toBe(!initial)
  })
})
