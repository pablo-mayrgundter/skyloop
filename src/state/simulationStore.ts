import {create} from 'zustand'
import {DEFAULT_SIM_SPEED} from '../config/simulation'
import {OverlayFlags, SimulationMode, SkyloopModuleId} from '../types'

type SimulationState = {
  time: number
  isPlaying: boolean
  simSpeed: number
  mode: SimulationMode
  selectedModule: SkyloopModuleId
} & OverlayFlags

type SimulationActions = {
  play: () => void
  pause: () => void
  reset: () => void
  setTime: (time: number) => void
  setSimSpeed: (speed: number) => void
  setMode: (mode: SimulationMode) => void
  setSelectedModule: (module: SkyloopModuleId) => void
  toggleOverlay: (key: keyof OverlayFlags) => void
}

const defaultOverlay: OverlayFlags = {
  showAtmosphere: true,
  showWind: false,
  showTethers: true,
  showSections: true,
  showGliders: true,
  showDayNight: true
}

export const useSimulationStore = create<SimulationState & SimulationActions>((set) => ({
  time: 0,
  isPlaying: false,
  simSpeed: DEFAULT_SIM_SPEED,
  mode: 'birdseye',
  selectedModule: 'main',
  ...defaultOverlay,
  play: () => set({ isPlaying: true }),
  pause: () => set({ isPlaying: false }),
  reset: () => set({ time: 0 }),
  setTime: (time) => set({ time }),
  setSimSpeed: (speed) => set({ simSpeed: speed }),
  setMode: (mode) => set({ mode }),
  setSelectedModule: (module) => set({ selectedModule: module }),
  toggleOverlay: (key) =>
    set((state) => ({
      [key]: !state[key]
    }))
}))
