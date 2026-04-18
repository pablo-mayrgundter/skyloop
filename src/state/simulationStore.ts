import { create } from 'zustand'
import { DEFAULT_SIM_SPEED, LOS_ANGELES, SAN_FRANCISCO } from '../config/simulation'
import { OverlayFlags, SimulationMode, SkyloopModuleId } from '../types/simulation'

export type SimulationState = OverlayFlags & {
  time: number
  isPlaying: boolean
  simSpeed: number
  mode: SimulationMode
  selectedModule: SkyloopModuleId
  lastUpdate: number
  route: {
    start: { lat: number; lon: number }
    end: { lat: number; lon: number }
  }
  play: () => void
  pause: () => void
  reset: () => void
  setTime: (time: number) => void
  setSimSpeed: (speed: number) => void
  setMode: (mode: SimulationMode) => void
  setSelectedModule: (module: SkyloopModuleId) => void
  toggleOverlay: (key: keyof OverlayFlags) => void
}

const overlayDefaults: OverlayFlags = {
  showAtmosphere: true,
  showWind: false,
  showTethers: true,
  showSections: true,
  showGliders: true,
  showDayNight: true
}

export const useSimulationStore = create<SimulationState>((set, get) => ({
  ...overlayDefaults,
  time: 0,
  isPlaying: false,
  simSpeed: DEFAULT_SIM_SPEED,
  mode: 'birdseye',
  selectedModule: 'main',
  lastUpdate: Date.now(),
  route: {
    start: LOS_ANGELES,
    end: SAN_FRANCISCO
  },
  play: () => set({ isPlaying: true, lastUpdate: Date.now() }),
  pause: () => set({ isPlaying: false }),
  reset: () => set({ time: 0 }),
  setTime: (time) => set({ time }),
  setSimSpeed: (speed) => set({ simSpeed: speed }),
  setMode: (mode) => set({ mode }),
  setSelectedModule: (module) => set({ selectedModule: module }),
  toggleOverlay: (key) => set({ [key]: !get()[key] } as Partial<SimulationState>)
}))
