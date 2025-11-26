import {describe, expect, it} from 'vitest'
import {useSimulationStore} from './simulationStore'

describe('simulation store', () => {
  it('toggles overlays', () => {
    const toggle = useSimulationStore.getState().toggleOverlay
    toggle('showAtmosphere')
    expect(useSimulationStore.getState().showAtmosphere).toBe(false)
  })

  it('sets mode and module', () => {
    const {setMode, setSelectedModule} = useSimulationStore.getState()
    setMode('firstPerson')
    setSelectedModule('economics')
    expect(useSimulationStore.getState().mode).toBe('firstPerson')
    expect(useSimulationStore.getState().selectedModule).toBe('economics')
  })
})
